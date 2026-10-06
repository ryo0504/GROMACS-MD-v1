"""Mechanical integration only: synthetic charges are NOT scientific RESP results."""
import json
from pathlib import Path
import pytest
from mdworkbench.config import load_config
from mdworkbench.utils.subprocess import detect
from mdworkbench.gromacs.runner import detect_gromacs, run_stage
from mdworkbench.parameterization.amber import parameterize

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def parameterized(tmp_path):
    pytest.importorskip("rdkit")
    pytest.importorskip("parmed")
    for executable in ("antechamber", "parmchk2", "tleap", "acpype"):
        if not detect(executable)["available"]:
            pytest.skip(f"{executable} is not installed")
    from mdworkbench.chemistry.build import build_molecule, write_sdf
    mol = build_molecule("CCO", 123)
    structure = tmp_path / "input.sdf"
    write_sdf(mol, structure)
    charges = [0.0] * 9
    oxygen = next(a for a in mol.GetAtoms() if a.GetSymbol() == "O")
    hydrogen = next(a for a in oxygen.GetNeighbors() if a.GetSymbol() == "H")
    charges[oxygen.GetIdx()] = -0.4
    charges[hydrogen.GetIdx()] = 0.4
    source = tmp_path / "SYNTHETIC_TEST_ONLY.charges"
    source.write_text(" ".join(map(str, charges)))
    cfg = load_config(ROOT / "systems/example/system.yaml")
    cfg["resp"].update(source_file=str(source), source_kind="charges", atom_order_confirmed=True,
                       purpose="synthetic test fixture; NOT fitted RESP; no scientific inference")
    cfg["structure"]["parameterization_file"] = str(structure)
    return cfg, parameterize(cfg, tmp_path / "parameters")


@pytest.mark.integration
def test_amber_to_gromacs_pipeline(parameterized):
    cfg, output = parameterized
    data = json.loads((output / "metadata.json").read_text())
    assert data["validation"]["atom_count"] == 9
    assert abs(data["validation"]["total_charge"]) < 1e-4
    assert len(data["commands"]) == 4
    assert all(c["status"] == "completed" for c in data["commands"])
    assert list(output.rglob("*_GMX.top"))
    assert list(output.rglob("*_GMX.gro"))


@pytest.mark.integration
def test_all_gromacs_stages_and_input_tampering(parameterized, tmp_path):
    if not detect_gromacs()["available"]:
        pytest.skip("GROMACS is not installed")
    import parmed as pmd
    from mdworkbench.gromacs.prepare import prepare_run
    from mdworkbench.provenance.metadata import hash_file
    cfg, output = parameterized
    top = next(output.rglob("*_GMX.top"))
    gro = next(output.rglob("*_GMX.gro"))
    structure = pmd.load_file(str(gro))
    structure.box = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]
    periodic = tmp_path / "periodic.gro"
    structure.save(str(periodic))
    templates = tmp_path / "mdp"
    templates.mkdir()
    import re
    for source in (ROOT / "mdp").glob("*.mdp"):
        text = re.sub(r"(?m)^nsteps\s*=.*$", "nsteps = 4", source.read_text())
        text = re.sub(r"(?m)^nstxout-compressed\s*=.*$", "nstxout-compressed = 1", text)
        (templates / source.name).write_text(text)
    run = prepare_run(cfg, periodic, top, templates, tmp_path / "run")
    for stage in ("em", "nvt", "npt", "production"):
        run_stage(run, stage, mdrun_args=["-ntmpi", "1"])
    data = json.loads((run / "metadata.json").read_text())
    assert all(value == "completed" for value in data["stages"].values())
    assert len(data["commands"]) == 8
    assert len(data["processed_hashes"]) == 4
    assert (run / "production/production.xtc").is_file()
    with pytest.raises(FileExistsError):
        run_stage(run, "production")
    target = run / "input/mdp/md.mdp"
    target.write_text(target.read_text() + "; changed\n")
    with pytest.raises(ValueError, match="Input changed"):
        run_stage(run, "production")
    assert not (run / ".running").exists()
