from pathlib import Path
from datetime import datetime, timezone
from importlib import metadata as packages
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
from mdworkbench.utils.subprocess import detect


def hash_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_info(directory):
    def get(*args):
        p = subprocess.run(["git", "-C", str(directory), *args], capture_output=True, text=True, check=False)
        return p.stdout.strip() if p.returncode == 0 else None
    if not detect("git")["available"]:
        return {"commit": None, "dirty": None}
    status = get("status", "--porcelain")
    return {"commit": get("rev-parse", "HEAD"), "dirty": bool(status) if status is not None else None,
            "status": status}


def collect_metadata(inputs=(), seed=None, config=None):
    from mdworkbench.gromacs.runner import detect_gromacs
    distributions = {d.metadata["Name"]: d.version for d in packages.distributions() if d.metadata.get("Name")}
    torch_info = {"version": distributions.get("torch"), "cuda_available": None,
                  "cuda_build": None, "mps_available": None}
    try:
        import torch
        torch_info.update(cuda_available=torch.cuda.is_available(), cuda_build=torch.version.cuda,
                          mps_available=torch.backends.mps.is_available())
    except ImportError:
        pass
    conda_meta = Path(os.environ.get("CONDA_PREFIX", sys.prefix)) / "conda-meta"
    amber = []
    for record in conda_meta.glob("ambertools-*.json"):
        data = json.loads(record.read_text())
        amber.append({k: data.get(k) for k in ("name", "version", "build", "subdir", "channel")})
    gpu = detect("nvidia-smi", ("--query-gpu=name,driver_version", "--format=csv"))
    if platform.system() == "Darwin":
        gpu = detect("system_profiler", ("SPDisplaysDataType", "-json"))
    return {"schema_version": 1, "timestamp": datetime.now(timezone.utc).isoformat(),
            "hostname": socket.gethostname(), "os": platform.platform(), "architecture": platform.machine(),
            "python": sys.version, "python_executable": sys.executable,
            "git": git_info(Path.cwd()),
            "micromamba": {"executable": detect("micromamba").get("executable"),
                "environment": os.environ.get("CONDA_DEFAULT_ENV"), "prefix": os.environ.get("CONDA_PREFIX"),
                "root_prefix": os.environ.get("MAMBA_ROOT_PREFIX")},
            "packages": distributions, "gromacs": detect_gromacs(), "ambertools": amber,
            "acpype": distributions.get("acpype"), "mace": distributions.get("mace-torch"),
            "pytorch": torch_info, "gpu": gpu,
            "cuda": detect("nvcc", ("--version",)), "random_seed": seed,
            "input_hashes": {str(Path(p).resolve()): hash_file(p) for p in inputs},
            "scientific_config": config, "commands": []}


def write_metadata(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str) + "\n")
    temporary.replace(path)
