#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
cd "$ROOT"
exec "$MM" run -n gmx-md jupyter lab "$ROOT/examples/ethanol_nvt_comparison/notebooks/ethanol_comparison.ipynb" "$@"
