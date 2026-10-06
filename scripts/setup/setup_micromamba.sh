#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
if [ ! -x "$MM" ]; then
  case "$(uname -s)/$(uname -m)" in
    Darwin/arm64) platform=osx-arm64;;
    Darwin/x86_64) platform=osx-64;;
    Linux/x86_64) platform=linux-64;;
    Linux/aarch64) platform=linux-aarch64;;
    *) echo 'Unsupported platform; install Micromamba explicitly.' >&2; exit 1;;
  esac
  mkdir -p "$ROOT/.local/bootstrap"
  curl -fL --retry 3 "https://micro.mamba.pm/api/micromamba/$platform/latest" -o "$ROOT/.local/bootstrap/micromamba.tar.bz2"
  tar -xjf "$ROOT/.local/bootstrap/micromamba.tar.bz2" -C "$ROOT/.local" bin/micromamba
fi
platform="$($MM info --json | /usr/bin/env python3 -c 'import json,sys; print(json.load(sys.stdin)["platform"])')"
if [ -d "$MAMBA_ROOT_PREFIX/envs/gmx-md" ]; then
  echo 'Existing gmx-md preserved. To change it, use micromamba install explicitly.'
elif [ -f "$ROOT/env/locks/$platform.lock" ]; then
  "$MM" create -y -n gmx-md -f "$ROOT/env/locks/$platform.lock"
else
  "$MM" create -y -n gmx-md --strict-channel-priority -f "$ROOT/env/environment.yml"
fi
"$MM" run -n gmx-md python -c 'import subprocess,sys,yaml; cfg=yaml.safe_load(open(sys.argv[1])); specs=[p for d in cfg["dependencies"] if isinstance(d,dict) for p in d.get("pip",[])]; subprocess.run([sys.executable,"-m","pip","install",*specs],check=True) if specs else None' "$ROOT/env/environment.yml"
"$MM" run -n gmx-md python -m pip install --no-deps --no-build-isolation -e "$ROOT"
"$MM" run -n gmx-md python "$ROOT/scripts/setup/register_kernel.py"
"$MM" env export -n gmx-md --explicit > "$ROOT/env/locks/$platform.lock"
echo 'Ready. See README for activation (no shell rc file was edited).'
