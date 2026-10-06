import argparse
import json
from pathlib import Path
import sys
from mdworkbench.config import load_config, resolve_path, positive_number


def build(cfg, output):
    from mdworkbench.chemistry.build import build_molecule, read_molecule, write_sdf
    from mdworkbench.chemistry.convert import to_ase, update_coordinates
    from mdworkbench.provenance.metadata import collect_metadata, write_metadata
    from ase.io import write
    structure = cfg["structure"]
    inputs = [Path(cfg["_config_path"])]
    if structure.get("source_file"):
        source = resolve_path(cfg, structure["source_file"])
        molecule = read_molecule(source)
        inputs.append(source)
    else:
        molecule = build_molecule(structure["smiles"], cfg["system"]["random_seed"])
    from rdkit import Chem
    if Chem.GetFormalCharge(molecule) != cfg["system"]["net_charge"]:
        raise ValueError("Structure formal charge does not match system.net_charge")
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    atoms = to_ase(molecule)
    mace = cfg.get("mace", {})
    if mace.get("use_geometry_optimization"):
        from mdworkbench.chemistry.mace import optimize
        if not mace.get("model_file"):
            raise ValueError("Explicit MACE model_file required")
        model = resolve_path(cfg, mace["model_file"])
        inputs.append(model)
        atoms = optimize(atoms, model, mace["device"], positive_number(mace["fmax_eV_per_A"], "fmax"),
                         int(positive_number(mace["steps"], "steps")), root / "optimization.traj")
        molecule = update_coordinates(molecule, atoms)
    write_sdf(molecule, root / "molecule.sdf")
    write(root / "molecule.xyz", atoms)
    metadata = collect_metadata(inputs, cfg["system"]["random_seed"], cfg)
    metadata["commands"] = [{"command": sys.argv, "status": "completed"}]
    write_metadata(root / "metadata.json", metadata)
    return root


def main():
    parser = argparse.ArgumentParser(description="GROMACS MD Workbench")
    sub = parser.add_subparsers(dest="operation", required=True)
    sub.add_parser("doctor")
    for name in ("build", "parameterize", "prepare-run"):
        p = sub.add_parser(name)
        p.add_argument("--config", required=True)
        p.add_argument("--output", required=True)
        if name == "prepare-run":
            p.add_argument("--structure", required=True)
            p.add_argument("--topology", required=True)
            p.add_argument("--mdp-dir", default="mdp")
            p.add_argument("--include-dir", action="append", default=[])
    p = sub.add_parser("run")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--stage", required=True, choices=("em", "nvt", "npt", "production"))
    p.add_argument("--platform", help="Optional host execution YAML; no implicit local config")
    args = parser.parse_args()
    try:
        if args.operation == "doctor":
            from mdworkbench.doctor import main as doctor
            doctor()
            return
        if args.operation == "run":
            from mdworkbench.gromacs.runner import run_stage
            import yaml
            host = yaml.safe_load(Path(args.platform).read_text()) if args.platform else {}
            if host.get("gromacs_executable"):
                import os
                os.environ["GMX_BIN"] = host["gromacs_executable"]
            for key in ("launcher", "mdrun_args"):
                if not isinstance(host.get(key, []), list) or not all(isinstance(x, str) for x in host.get(key, [])):
                    raise ValueError(f"platform.{key} must be a list of strings")
            run_stage(args.run_dir, args.stage, host.get("launcher", []), host.get("mdrun_args", []))
            return
        cfg = load_config(args.config)
        if args.operation == "build":
            result = build(cfg, args.output)
        elif args.operation == "parameterize":
            from mdworkbench.parameterization.amber import parameterize
            result = parameterize(cfg, args.output)
        else:
            from mdworkbench.gromacs.prepare import prepare_run
            result = prepare_run(cfg, args.structure, args.topology, args.mdp_dir, args.output, args.include_dir)
        print(result)
    except (ValueError, RuntimeError, OSError, KeyError) as error:
        parser.exit(2, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
