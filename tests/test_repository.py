from pathlib import Path
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_output_ignored_and_scientific_inputs_allowed():
    from mdworkbench.utils.subprocess import detect
    if not detect("git")["available"]:
        pytest.skip("git unavailable")
    ignored = ["runs/result.gro", "result.xtc", "result.trr", "result.edr", "result.cpt", "result.tpr",
               ".local/model.model", ".micromamba/envs/gmx-md/bin/python", ".ipynb_checkpoints/example.ipynb"]
    allowed = ["mdp/em.mdp", "systems/example/parameters/final.mol2", "systems/example/parameters/system.top",
               "systems/example/parameters/molecule.itp", "systems/example/structures/system.gro", "systems/example/parameters/resp.in"]
    for files, expected in ((ignored, 0), (allowed, 1)):
        for file in files:
            p = subprocess.run(["git", "check-ignore", "--no-index", file], cwd=ROOT, capture_output=True)
            # Downloaded source archives may lack .git; do not misclassify Git errors as successful validation.
            if p.returncode == 128:
                pytest.skip("repository not initialized")
            assert p.returncode == expected, file


def test_notebook_structure():
    nbformat = pytest.importorskip("nbformat")
    notebooks = list((ROOT / "notebooks").glob("*.ipynb"))
    assert len(notebooks) == 4
    for path in notebooks:
        notebook = nbformat.read(path, as_version=4)
        nbformat.validate(notebook)
        assert notebook.metadata.kernelspec.name == "gmx-md"
