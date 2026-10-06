#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
if [ -x "$MM" ] && [ -d "$MAMBA_ROOT_PREFIX/envs/gmx-md" ]; then
  exec "$MM" run -n gmx-md python -m mdworkbench.doctor "$@"
elif command -v python3 >/dev/null 2>&1; then
  export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
  exec python3 -m mdworkbench.doctor "$@"
else
  echo 'Python未導入: run scripts/setup/setup_micromamba.sh'; exit 0
fi
