#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../setup/common.sh"
exec python -m mdworkbench.cli run --stage nvt "$@"
