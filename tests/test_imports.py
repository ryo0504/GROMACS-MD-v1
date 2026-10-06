import importlib
import pytest


def test_core_imports():
    for module in ("mdworkbench", "mdworkbench.config", "mdworkbench.gromacs.runner", "mdworkbench.parameterization.resp"):
        importlib.import_module(module)


@pytest.mark.parametrize("module", ["rdkit", "ase", "MDAnalysis", "mdtraj", "parmed", "mace.calculators", "torch", "jupyterlab"])
def test_scientific_imports(module):
    pytest.importorskip(module)
