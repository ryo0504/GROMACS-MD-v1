import os
from pathlib import Path
import pytest
from mdworkbench.gromacs.runner import detect_gromacs
from mdworkbench.utils.subprocess import run


def test_missing_gromacs(monkeypatch):
    monkeypatch.setenv("PATH", "")
    monkeypatch.delenv("GMX_BIN", raising=False)
    assert not detect_gromacs()["available"]


def test_explicit_executable_with_spaces(tmp_path, monkeypatch):
    script = tmp_path / "fake gmx"
    script.write_text("#!/bin/sh\necho 'GROMACS version: test'\n")
    script.chmod(0o755)
    monkeypatch.setenv("GMX_BIN", str(script))
    assert detect_gromacs()["usable"]


def test_bad_override_does_not_fall_back(monkeypatch):
    monkeypatch.setenv("GMX_BIN", "/missing/gmx")
    assert not detect_gromacs()["usable"]


@pytest.mark.integration
def test_gromacs_version():
    gmx = detect_gromacs()
    if not gmx["available"]: pytest.skip("External GROMACS missing")
    assert gmx["usable"]
    assert "GROMACS version" in gmx["output"]


def test_failed_external_command(tmp_path):
    with pytest.raises(RuntimeError, match="Command failed"):
        run(["/bin/sh", "-c", "exit 7"], tmp_path, tmp_path / "command.log")
