#!/usr/bin/env bash
set -u
source "$(dirname "$0")/common.sh"
if [ -n "${GMX_BIN:-}" ]; then
  if command -v "$GMX_BIN" >/dev/null 2>&1; then "$GMX_BIN" --version; else echo 'GMX_BIN is unavailable'; fi
else
  found=false
  for exe in gmx gmx_mpi; do
    if command -v "$exe" >/dev/null 2>&1; then command -v "$exe"; "$exe" --version; found=true; fi
  done
  if [ "$found" = false ]; then echo 'GROMACS未設定: see env/README.md'; fi
fi
