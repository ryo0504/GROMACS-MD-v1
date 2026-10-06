import importlib
import json
import os
import platform
import sys
from mdworkbench.gromacs.runner import detect_gromacs
from mdworkbench.utils.subprocess import detect


def report():
    modules = {}
    for name in ("numpy", "scipy", "pandas", "matplotlib", "rdkit", "ase", "MDAnalysis", "mdtraj", "parmed",
                 "jupyterlab", "ipykernel", "yaml", "pytest", "mace", "torch", "acpype"):
        try:
            imported = importlib.import_module(name)
            modules[name] = {"available": True, "version": getattr(imported, "__version__", "unknown")}
        except Exception as error:
            modules[name] = {"available": False, "error": str(error)}
    device = {}
    if modules["torch"]["available"]:
        import torch
        device = {"cuda_available": torch.cuda.is_available(), "cuda_build": torch.version.cuda,
                  "mps_available": torch.backends.mps.is_available()}
    return {"os": platform.platform(), "architecture": platform.machine(), "python": sys.version,
            "prefix": sys.prefix, "environment": os.environ.get("CONDA_DEFAULT_ENV"),
            "micromamba": detect("micromamba", ("--version",)), "modules": modules,
            "ambertools": {name: detect(name) for name in ("antechamber", "parmchk2", "tleap", "resp", "espgen")}, "acpype": detect("acpype"),
            "gromacs": detect_gromacs(), "torch_device": device}


def main():
    data = report()
    if "--json" in sys.argv:
        print(json.dumps(data, indent=2))
    else:
        print(f"[INFO] {data['os']} / {data['architecture']} / Python {sys.version.split()[0]}")
        for name, entry in {"micromamba": data["micromamba"], **data["modules"],
                            **data["ambertools"], "ACPYPE CLI": data["acpype"]}.items():
            print(f"[{'OK' if entry['available'] else 'MISSING'}] {name}")
            if not entry["available"]:
                print("  Fix: scripts/setup/setup_micromamba.sh; activate gmx-md; MACE: scripts/setup/install_mace.sh")
                if "error" in entry:
                    print(" ", entry["error"])
        print(f"[{'OK' if data['gromacs'].get('usable') else 'INFO'}] GROMACS: {data['gromacs']['output']}")
        print("[INFO] torch device:", data["torch_device"])
    # Inventory is advisory: missing GROMACS never invalidates the Python environment.


if __name__ == "__main__":
    main()
