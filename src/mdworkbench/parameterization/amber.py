from pathlib import Path
import shutil
from mdworkbench.chemistry.build import read_molecule
from mdworkbench.config import resolve_path
from mdworkbench.parameterization.resp import require_source, fit_or_load
from mdworkbench.parameterization.validation import validate_frcmod, topology_summary
from mdworkbench.provenance.metadata import collect_metadata, write_metadata
from mdworkbench.utils.subprocess import detect, run


def detect_ambertools():
    return {name: detect(name) for name in ("antechamber", "parmchk2", "tleap", "resp", "espgen")}


def parameterize(cfg, output):
    source = require_source(cfg)  # fail before dependency checks; never choose another charge model
    structure_path = resolve_path(cfg, cfg["structure"]["parameterization_file"])
    molecule = read_molecule(structure_path)
    n = molecule.GetNumAtoms()
    missing = [k for k, v in detect_ambertools().items() if not v["available"] and k in {"antechamber", "parmchk2", "tleap"}]
    if missing or not detect("acpype")["available"]:
        raise RuntimeError(f"Missing AmberTools/ACPYPE: {missing}")
    work = Path(output).resolve()
    work.mkdir(parents=True, exist_ok=False)
    # LEaP reads simple local filenames even when the repository path contains spaces.
    input_structure = work / ("input" + structure_path.suffix.lower())
    shutil.copy2(structure_path, input_structure)
    inputs = [Path(cfg["_config_path"]), structure_path, source]
    inputs += [resolve_path(cfg, p) for p in cfg["resp"].get("fit_inputs", [])]
    metadata = collect_metadata(inputs, cfg["system"]["random_seed"], cfg)
    write_metadata(work / "metadata.json", metadata)
    commands = metadata["commands"]

    def execute(command, log):
        record = {"command": [str(x) for x in command], "status": "started"}
        commands.append(record)
        write_metadata(work / "metadata.json", metadata)
        try:
            run(command, work, work / log)
            record["status"] = "completed"
        except Exception:
            record["status"] = "failed"
            raise
        finally:
            write_metadata(work / "metadata.json", metadata)

    charges = fit_or_load(cfg, work, n, execute)
    ff = cfg["force_field"]["family"]
    execute(["antechamber", "-i", input_structure.name, "-fi", "sdf" if input_structure.suffix == ".sdf" else "mol2",
             "-o", "final.mol2", "-fo", "mol2", "-at", ff, "-c", "rc", "-cf", charges.name,
             "-nc", str(cfg["system"]["net_charge"]), "-m", str(cfg["system"]["multiplicity"]), "-s", "2"], "antechamber.log")
    execute(["parmchk2", "-i", "final.mol2", "-f", "mol2", "-o", "parameters.frcmod", "-s", ff], "parmchk2.log")
    validate_frcmod(work / "parameters.frcmod")
    (work / "tleap.in").write_text(f"source leaprc.{ff}\nloadamberparams parameters.frcmod\nMOL = loadmol2 final.mol2\ncheck MOL\nsaveamberparm MOL system.prmtop system.inpcrd\nquit\n")
    execute(["tleap", "-f", "tleap.in"], "tleap.log")
    metadata["validation"] = topology_summary(work / "system.prmtop", work / "system.inpcrd",
                                              cfg["system"]["net_charge"], n)
    execute(["acpype", "-p", "system.prmtop", "-x", "system.inpcrd"], "acpype.log")
    metadata["output_hashes"] = collect_metadata(list(work.rglob("*.mol2")) + list(work.rglob("*.frcmod")) +
         list(work.rglob("*.prmtop")) + list(work.rglob("*.top")) + list(work.rglob("*.itp")) + list(work.rglob("*.gro")))["input_hashes"]
    write_metadata(work / "metadata.json", metadata)
    return work
