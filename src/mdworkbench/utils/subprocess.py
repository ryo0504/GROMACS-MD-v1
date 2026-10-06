"""Argument-vector execution only; commands and exit status remain reviewable."""
import os
from pathlib import Path
import shutil
import subprocess


def detect(name, args=(), timeout=20):
    executable = shutil.which(str(name))
    if not executable:
        return {"available": False, "executable": None, "output": "未導入"}
    result = {"available": True, "executable": executable}
    if args:
        try:
            p = subprocess.run([executable, *args], capture_output=True, text=True,
                               timeout=timeout, check=False)
            result.update(output=p.stdout + p.stderr, returncode=p.returncode)
        except (OSError, subprocess.TimeoutExpired) as error:
            result.update(output=str(error), returncode=None)
    return result


def run(command, cwd, log=None, timeout=None):
    command = [os.fspath(arg) for arg in command]
    if log:
        with Path(log).open("w") as handle:
            p = subprocess.run(command, cwd=cwd, stdout=handle,
                               stderr=subprocess.STDOUT, timeout=timeout, check=False)
    else:
        p = subprocess.run(command, cwd=cwd, timeout=timeout, check=False)
    if p.returncode:
        raise RuntimeError(f"Command failed ({p.returncode}): {command}; log={log}")
    return p
