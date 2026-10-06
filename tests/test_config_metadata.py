from pathlib import Path
import json
import pytest
from mdworkbench.config import load_config, resolve_path
from mdworkbench.provenance.metadata import collect_metadata, hash_file, write_metadata

ROOT = Path(__file__).resolve().parents[1]


def test_config_loading():
    cfg = load_config(ROOT / "systems/example/system.yaml")
    assert cfg["force_field"]["family"] == "gaff2"
    assert resolve_path(cfg, "structures/molecule.sdf") == ROOT / "systems/example/structures/molecule.sdf"


@pytest.mark.parametrize("mutation", ["charge", "seed", "basis", "force_field"])
def test_invalid_config(tmp_path, mutation):
    import yaml
    cfg = load_config(ROOT / "systems/example/system.yaml")
    if mutation == "charge": cfg["system"]["net_charge"] = "0"
    elif mutation == "seed": cfg["system"]["random_seed"] = -1
    elif mutation == "basis": cfg["resp"]["basis"] = None
    else: cfg["force_field"]["charge_model"] = "bcc"
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError): load_config(path)


def test_metadata_hash_and_roundtrip(tmp_path):
    file = tmp_path / "input.mdp"
    file.write_text("gen-seed = 42\n")
    data = collect_metadata([file], 42)
    assert data["input_hashes"][str(file)] == hash_file(file)
    assert data["random_seed"] == 42
    assert {"git", "hostname", "architecture", "gpu", "gromacs", "ambertools", "pytorch", "micromamba"} <= data.keys()
    output = tmp_path / "metadata.json"
    write_metadata(output, data)
    assert json.loads(output.read_text()) == data
    assert not output.with_suffix(".json.tmp").exists()
    file.write_text("changed")
    assert hash_file(file) != data["input_hashes"][str(file)]
