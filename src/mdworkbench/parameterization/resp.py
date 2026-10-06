"""External RESP ingress. Fitting constraints and QM settings are user supplied."""
from pathlib import Path
import math
from mdworkbench.config import resolve_path
from mdworkbench.utils.subprocess import run, detect


def require_source(cfg):
    resp = cfg["resp"]
    if not resp.get("source_file"):
        raise ValueError("RESP source_file is required. AM1-BCC fallback is prohibited.")
    source = resolve_path(cfg, resp["source_file"])
    if not source.is_file():
        raise ValueError(f"RESP source does not exist: {source}")
    if resp.get("atom_order_confirmed") is not True:
        raise ValueError("Set resp.atom_order_confirmed only after checking structure/QM/charge atom order")
    if resp.get("source_kind") not in {"charges", "esp", "gaussian_output"}:
        raise ValueError("source_kind must be charges, esp or gaussian_output")
    return source


def read_charges(path, atom_count, net_charge, tolerance=1e-4):
    values = [float(token.replace("D", "E")) for token in Path(path).read_text().split()]
    if len(values) != atom_count or not all(math.isfinite(q) for q in values):
        raise ValueError("RESP charge count/finite-value check failed")
    if abs(sum(values) - net_charge) > tolerance:
        raise ValueError(f"RESP net charge mismatch: {sum(values)} vs {net_charge}")
    return values


def fit_or_load(cfg, work, atom_count, execute=None):
    source = require_source(cfg)
    resp = cfg["resp"]
    if execute is None:
        execute = lambda command, log: run(command, work, work / log)
    if resp["source_kind"] == "charges":
        result = source
    else:
        stages = resp.get("fit_inputs")
        if not stages:
            raise ValueError("RESP fit_inputs must contain user-reviewed restraint/equivalence controls")
        esp = source
        if resp["source_kind"] == "gaussian_output":
            if not detect("espgen")["available"]:
                raise RuntimeError("espgen is missing (AmberTools)")
            esp = work / "qm.esp"
            execute(["espgen", "-i", source, "-o", esp], "espgen.log")
        if not detect("resp")["available"]:
            raise RuntimeError("resp executable is missing (AmberTools)")
        previous = None
        for i, fit_input in enumerate(stages, 1):
            fit_input = resolve_path(cfg, fit_input)
            fit_input.resolve(strict=True)
            result = work / f"resp_stage{i}.qout"
            command = ["resp", "-O", "-i", fit_input, "-o", work / f"resp_stage{i}.out",
                       "-e", esp, "-t", result, "-p", work / f"resp_stage{i}.pch"]
            if previous:
                command += ["-q", previous]
            execute(command, f"resp_stage{i}.log")
            previous = result
    values = read_charges(result, atom_count, cfg["system"]["net_charge"])
    target = work / "resp.charges"
    target.write_text("\n".join(f"{q:.10f}" for q in values) + "\n")
    return target
