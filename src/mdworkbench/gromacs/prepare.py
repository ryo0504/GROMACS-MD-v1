from pathlib import Path
import os
import re
import shutil
import tempfile
from mdworkbench.config import positive_number
from mdworkbench.gromacs.validation import INCLUDE, validate_gro, mdp_values
from mdworkbench.provenance.metadata import collect_metadata, hash_file, write_metadata


def snapshot_topology(topology, destination, search_paths=()):
    """Flatten and rewrite the complete include graph (including conditional includes)."""
    destination.mkdir()
    visited = {}

    def copy(source, filename=None):
        source = Path(source).resolve(strict=True)
        if source in visited:
            return visited[source]
        filename = filename or f"{hash_file(source)[:12]}_{source.name}"
        visited[source] = filename
        text = source.read_text()

        def include(match):
            name = match.group(1)
            candidates = [source.parent / name, *[Path(p) / name for p in search_paths]]
            resolved = next((p for p in candidates if p.is_file()), None)
            if resolved is None:
                raise ValueError(f"Missing topology include {name} referenced by {source}; supply --include-dir")
            return '#include "' + copy(resolved) + '"'

        (destination / filename).write_text(INCLUDE.sub(include, text))
        return filename

    copy(topology, "system.top")
    return list(visited)


def prepare_run(cfg, structure, topology, mdp_dir, output, include_dirs=()):
    structure, topology = Path(structure).resolve(strict=True), Path(topology).resolve(strict=True)
    count = validate_gro(structure)
    # Resolve caller-provided include dirs before ParmEd parsing; validate in a temporary snapshot.
    from mdworkbench.parameterization.validation import topology_summary
    import parmed as pmd
    from parmed.gromacs import GromacsGroFile
    root = Path(output).resolve()
    with tempfile.TemporaryDirectory(prefix="mdworkbench-topology-") as temporary:
        staged = Path(temporary) / "topology"
        sources = snapshot_topology(topology, staged, include_dirs)
        top = pmd.load_file(str(staged / "system.top"))
        gro = GromacsGroFile.parse(str(structure))
        if len(top.atoms) != count or len(top.residues) != len(gro.residues):
            raise ValueError("Coordinate/topology atom or residue count mismatch")
        if [a.name for a in top.atoms] != [a.name for a in gro.atoms]:
            raise ValueError("Coordinate/topology atom order mismatch")
        summary = topology_summary(staged / "system.top", net_charge=cfg["system"]["net_charge"], atom_count=count)
        root.mkdir(parents=True, exist_ok=False)
        (root / "input/mdp").mkdir(parents=True)
        shutil.copytree(staged, root / "input/topology")
    shutil.copy2(structure, root / "input/structure.gro")
    shutil.copy2(cfg["_config_path"], root / "input/system.yaml")
    temperature = positive_number(cfg["gromacs"]["temperature_K"], "temperature_K")
    pressure = positive_number(cfg["gromacs"]["pressure_bar"], "pressure_bar")
    for name in ("em", "nvt", "npt", "md"):
        source = Path(mdp_dir) / f"{name}.mdp"
        text = source.read_text()
        replacements = {"ref-t": temperature, "gen-temp": temperature, "ref-p": pressure,
                        "gen-seed": cfg["system"]["random_seed"]}
        for key, value in replacements.items():
            text = re.sub(rf"(?m)^({key}\s*=).*?$", rf"\g<1> {value}", text)
        target = root / "input/mdp" / source.name
        target.write_text(text)
        if name == "nvt" and mdp_values(target).get("gen-seed") != str(cfg["system"]["random_seed"]):
            raise ValueError("NVT MDP must expose explicit gen-seed")
    inputs = [structure, Path(cfg["_config_path"]), *sources, *Path(mdp_dir).glob("*.mdp")]
    metadata = collect_metadata(inputs, cfg["system"]["random_seed"], cfg)
    metadata["validation"] = summary
    metadata["snapshot_hashes"] = {str(p.relative_to(root)): hash_file(p) for p in (root / "input").rglob("*") if p.is_file()}
    write_metadata(root / "metadata.json", metadata)
    return root
