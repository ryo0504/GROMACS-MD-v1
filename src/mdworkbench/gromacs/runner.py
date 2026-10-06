from pathlib import Path
import json
import os
import shutil
from datetime import datetime, timezone
from mdworkbench.gromacs.validation import validate_gro
from mdworkbench.provenance.metadata import hash_file, write_metadata, collect_metadata
from mdworkbench.utils.subprocess import detect, run


def detect_gromacs():
    override = os.environ.get("GMX_BIN")
    names = [override] if override else ["gmx", "gmx_mpi"]
    for name in names:
        result = detect(name, ("--version",))
        if result["available"]:
            result["usable"] = result.get("returncode") == 0
            return result
    return {"available": False, "usable": False, "executable": None, "output": "GROMACS未設定"}


def run_stage(run_dir, stage, launcher=(), mdrun_args=()):
    root = Path(run_dir).resolve()
    order = ["em", "nvt", "npt", "production"]
    if stage not in order:
        raise ValueError("Unknown stage")
    gmx = detect_gromacs()
    if not gmx.get("usable"):
        raise RuntimeError(f"GROMACS is unavailable/unusable: {gmx['output']}")
    metadata_path = root / "metadata.json"
    # Atomic mkdir is a per-run execution lock; never run two stages in parallel.
    lock = root / ".running"
    lock.mkdir()
    try:
        metadata = json.loads(metadata_path.read_text())
        for relative, digest in metadata["snapshot_hashes"].items():
            if hash_file(root / relative) != digest:
                raise ValueError(f"Input changed since prepare-run: {relative}")
        i = order.index(stage)
        previous = order[i - 1] if i else None
        coordinates = root / previous / f"{previous}.gro" if previous else root / "input/structure.gro"
        if previous and metadata.get("stages", {}).get(previous) != "completed":
            raise ValueError(f"Stage {previous} must complete first")
        validate_gro(coordinates)
        work = root / stage
        work.mkdir(exist_ok=False)  # prevents accidental checkpoint/output overwrite
        mdp = root / "input/mdp" / ("md.mdp" if stage == "production" else f"{stage}.mdp")
        command = [gmx["executable"], "grompp", "-f", mdp, "-c", coordinates,
                   "-p", root / "input/topology/system.top", "-o", work / f"{stage}.tpr",
                   "-po", work / "processed.mdp", "-pp", work / "processed.top"]
        if stage in {"npt", "production"}:
            checkpoint = root / previous / f"{previous}.cpt"
            checkpoint.resolve(strict=True)
            command += ["-t", checkpoint]
        metadata["gromacs"] = gmx
        metadata.setdefault("stage_environments", {})[stage] = collect_metadata()

        def execute(argv, log):
            entry = {"command": [str(a) for a in argv], "cwd": str(work),
                     "timestamp": datetime.now(timezone.utc).isoformat(), "status": "started"}
            metadata["commands"].append(entry)
            write_metadata(metadata_path, metadata)
            try:
                run(argv, work, work / log)
                entry["status"] = "completed"
            except Exception:
                entry["status"] = "failed"
                metadata.setdefault("stages", {})[stage] = "failed"
                raise
            finally:
                write_metadata(metadata_path, metadata)

        execute(command, "grompp.log")  # no -maxwarn: warnings remain blockers
        metadata.setdefault("processed_hashes", {})[stage] = {
            name: hash_file(work / name) for name in ("processed.mdp", "processed.top", f"{stage}.tpr")}
        execute([*launcher, gmx["executable"], "mdrun", "-deffnm", stage, *mdrun_args], "mdrun.log")
        validate_gro(work / f"{stage}.gro")
        metadata.setdefault("stages", {})[stage] = "completed"
        write_metadata(metadata_path, metadata)
    finally:
        lock.rmdir()
