"""Shared, explicit file/process handling for the server preparation bundle."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BOHR_A = 0.529177210903


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def execute(args, cwd, log, stdin=None):
    """Never use a shell, never hide failure, never silently overwrite a log."""
    args = [str(x) for x in args]
    if not shutil.which(args[0]) and not Path(args[0]).is_file():
        raise RuntimeError(f"Executable unavailable: {args[0]}")
    log = Path(log)
    with log.open("x") as handle:
        result = subprocess.run(args, cwd=cwd, input=stdin, text=True,
                                stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {args}; see {log}")
    return {"argv": args, "returncode": result.returncode, "log": str(log)}


def read_esp(path):
    """Amber RESP fixed-width ESP, coordinates in bohr, potential in a.u."""
    import numpy as np
    lines = Path(path).read_text().splitlines()
    n, p = int(lines[0][:5]), int(lines[0][5:10])
    if n <= 0 or p <= 0 or len(lines) != 1 + n + p:
        raise ValueError(f"Invalid/incomplete ESP dimensions: {path}")
    xyz = np.array([[float(l[17+16*i:33+16*i].replace('D', 'E'))
                     for i in range(3)] for l in lines[1:1+n]])
    data = np.array([[float(l[1+16*i:17+16*i].replace('D', 'E'))
                      for i in range(4)] for l in lines[1+n:]])
    if not np.isfinite(xyz).all() or not np.isfinite(data).all():
        raise ValueError("Nonfinite ESP data")
    if np.min(np.linalg.norm(data[:, None, 1:] - xyz[None], axis=2)) < .2:
        raise ValueError("ESP grid points dangerously close to nuclei")
    return xyz, data[:, 1:], data[:, 0]


def esp_error(path, charges):
    import numpy as np
    xyz, grid, observed = read_esp(path)
    predicted = np.zeros(len(grid))
    for start in range(0, len(grid), 2048):
        d = np.linalg.norm(grid[start:start+2048, None] - xyz[None], axis=2)
        predicted[start:start+2048] = (1 / d) @ np.asarray(charges)
    error = predicted - observed
    return {"n_points": len(grid), "rmse_au": float(np.sqrt(np.mean(error**2))),
            "max_abs_error_au": float(np.max(np.abs(error)))}
