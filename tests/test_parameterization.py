from pathlib import Path
import pytest
from mdworkbench.config import load_config
from mdworkbench.parameterization.resp import require_source, read_charges, fit_or_load
from mdworkbench.parameterization.amber import detect_ambertools, parameterize
from mdworkbench.parameterization.validation import validate_frcmod
from mdworkbench.utils.subprocess import detect

ROOT = Path(__file__).resolve().parents[1]


def test_resp_missing_never_falls_back(tmp_path):
    cfg = load_config(ROOT / "systems/example/system.yaml")
    with pytest.raises(ValueError, match="fallback"):
        parameterize(cfg, tmp_path / "parameters")
    assert not (tmp_path / "parameters").exists()


def test_resp_atom_order_confirmation(tmp_path):
    cfg = load_config(ROOT / "systems/example/system.yaml")
    source = tmp_path / "charges"
    source.write_text("0\n")
    cfg["resp"]["source_file"] = str(source)
    with pytest.raises(ValueError, match="atom_order"):
        require_source(cfg)


@pytest.mark.parametrize("content,count,charge", [("0.1 -0.1", 2, 0), ("0", 2, 0), ("nan 0", 2, 0), ("0.5 0.5", 2, 0)])
def test_resp_charge_validation(tmp_path, content, count, charge):
    file = tmp_path / "charges"
    file.write_text(content)
    if content == "0.1 -0.1":
        assert read_charges(file, count, charge) == [0.1, -0.1]
    else:
        with pytest.raises(ValueError): read_charges(file, count, charge)


def test_ambertools_detection_contract(monkeypatch):
    monkeypatch.setenv("PATH", "")
    assert all(not v["available"] for v in detect_ambertools().values())
    assert not detect("acpype")["available"]


@pytest.mark.integration
def test_installed_ambertools_and_acpype():
    if not all(v["available"] for v in detect_ambertools().values()) or not detect("acpype")["available"]:
        pytest.skip("AmberTools/ACPYPE not installed")
    assert detect("antechamber", ("-L",))["returncode"] == 0
    assert "acpype" in detect("acpype", ("-h",))["output"].lower()


def test_missing_parameters_rejected(tmp_path):
    file = tmp_path / "params.frcmod"
    file.write_text("ATTN, need revision")
    with pytest.raises(ValueError): validate_frcmod(file)
