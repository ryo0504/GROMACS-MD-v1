#!/usr/bin/env bash
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
if [ -f "$ROOT/.local/host.env" ]; then
  source "$ROOT/.local/host.env"
fi
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$ROOT/.micromamba}"
if [ -n "${MDWORKBENCH_MM:-}" ]; then
  MM="$MDWORKBENCH_MM"
elif command -v micromamba >/dev/null 2>&1; then
  MM="$(command -v micromamba)"
else
  MM="$ROOT/.local/bin/micromamba"
fi
export PATH="${MDWORKBENCH_HOST_BIN:+$MDWORKBENCH_HOST_BIN:}$ROOT/.local/gromacs/bin:$ROOT/.local/bin:$PATH"
