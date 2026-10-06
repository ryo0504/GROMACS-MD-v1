from pathlib import Path
import math
import re
import yaml


def load_config(filename):
    path = Path(filename).resolve()
    with path.open() as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise ValueError("Configuration must be a YAML mapping")
    system = cfg.get("system", {})
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", str(system.get("name", ""))):
        raise ValueError("system.name must be a safe identifier")
    for key in ("net_charge", "multiplicity", "random_seed"):
        value = system.get(key)
        if type(value) is not int:
            raise ValueError(f"system.{key} must be an explicit integer")
    if system["multiplicity"] < 1 or not 0 <= system["random_seed"] < 2**31:
        raise ValueError("Invalid multiplicity or reproducible random seed")
    ff = cfg.get("force_field", {})
    if ff.get("family") not in {"gaff", "gaff2"} or ff.get("charge_model") != "resp":
        raise ValueError("v1 requires explicit gaff/gaff2 and resp; no fallback")
    resp = cfg.get("resp", {})
    for key in ("qm_backend", "method", "basis"):
        if not isinstance(resp.get(key), str) or not resp[key].strip():
            raise ValueError(f"resp.{key} must be explicit")
    cfg["_config_path"] = str(path)
    return cfg


def resolve_path(cfg, value):
    path = Path(value)
    return (Path(cfg["_config_path"]).parent / path).resolve() if not path.is_absolute() else path.resolve()


def positive_number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value
