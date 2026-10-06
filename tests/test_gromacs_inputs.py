from pathlib import Path
import pytest
from mdworkbench.gromacs.prepare import snapshot_topology
from mdworkbench.gromacs.validation import validate_gro, mdp_values


def test_include_snapshot_and_missing_include(tmp_path):
    (tmp_path / "system.top").write_text('#include "molecule.itp"\n[ molecules ]\nMOL 1\n')
    (tmp_path / "molecule.itp").write_text("[ moleculetype ]\nMOL 3\n")
    sources = snapshot_topology(tmp_path / "system.top", tmp_path / "snapshot")
    assert len(sources) == 2
    assert "molecule.itp" in (tmp_path / "snapshot/system.top").read_text()
    (tmp_path / "molecule.itp").unlink()
    with pytest.raises(ValueError, match="Missing topology include"):
        snapshot_topology(tmp_path / "system.top", tmp_path / "bad")


def test_mdp_no_duplicates(tmp_path):
    file = tmp_path / "bad.mdp"
    file.write_text("gen_seed = 1\ngen-seed = 2\n")
    with pytest.raises(ValueError, match="Duplicate"): mdp_values(file)


def test_templates_are_explicit():
    root = Path(__file__).resolve().parents[1]
    assert mdp_values(root / "mdp/nvt.mdp")["gen-seed"] == "20261004"
    for stage in ("npt", "md"):
        values = mdp_values(root / f"mdp/{stage}.mdp")
        assert values["gen-vel"] == "no"
        assert values["continuation"] == "yes"


def test_empty_box_rejected(tmp_path):
    file = tmp_path / "bad.gro"
    file.write_text("test\n1\n    1MOL      C    1   0.000   0.000   0.000\n0 0 0\n")
    with pytest.raises(ValueError, match="periodic box"): validate_gro(file)
