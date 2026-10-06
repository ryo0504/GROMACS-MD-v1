from pathlib import Path


def evaluate(atoms, model, device="cpu"):
    """Use only an explicitly supplied model; never download/select a scientific model."""
    if device == "mps":
        raise ValueError("This float64 MACE workflow supports CPU/CUDA; MPS requires a separately reviewed float32 workflow")
    from mace.calculators import MACECalculator
    model = Path(model).resolve(strict=True)
    atoms = atoms.copy()
    atoms.calc = MACECalculator(model_paths=str(model), device=device, default_dtype="float64")
    return atoms, {"energy_eV": float(atoms.get_potential_energy()),
                   "forces_eV_per_A": atoms.get_forces().tolist()}


def optimize(atoms, model, device, fmax, steps, trajectory):
    from ase.optimize import BFGS
    atoms, _ = evaluate(atoms, model, device)
    converged = BFGS(atoms, trajectory=str(trajectory)).run(fmax=fmax, steps=steps)
    if not converged:
        raise RuntimeError("MACE geometry optimization did not converge")
    return atoms
