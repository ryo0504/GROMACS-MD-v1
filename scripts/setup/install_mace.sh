#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
echo "OS=$(uname -s) architecture=$(uname -m)"
if command -v nvcc >/dev/null 2>&1; then nvcc --version; else echo 'nvcc: 未導入'; fi
gpu=false
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
  gpu=true
fi
if ! "$MM" run -n gmx-md python -c 'import torch; print(torch.__version__, "CUDA", torch.version.cuda, "available", torch.cuda.is_available())'; then
  case "$(uname -s)" in
    Darwin) "$MM" install -y -n gmx-md -c conda-forge --strict-channel-priority 'pytorch>=2.2,<2.9';;
    Linux)
      if [ "$gpu" = true ]; then
        if [ -z "${TORCH_INDEX_URL:-}" ] || [ -z "${TORCH_VERSION:-}" ]; then
          echo 'Select a driver-compatible wheel using https://pytorch.org/get-started/locally/.' >&2
          echo 'Then set TORCH_INDEX_URL and TORCH_VERSION explicitly and rerun.' >&2
          exit 2
        fi
        "$MM" run -n gmx-md python -m pip install "torch==$TORCH_VERSION" --index-url "$TORCH_INDEX_URL"
      else
        "$MM" run -n gmx-md python -m pip install 'torch>=2.2,<2.9' --index-url https://download.pytorch.org/whl/cpu
      fi;;
    *) echo 'Unsupported PyTorch host'; exit 2;;
  esac
fi
if [ "$gpu" = true ]; then
  "$MM" run -n gmx-md python -c 'import torch; assert torch.cuda.is_available(), "Installed torch cannot use this NVIDIA host; select a compatible backend explicitly"'
fi
if [ -n "${MACE_PIP_LOCK:-}" ]; then
  "$MM" run -n gmx-md python -m pip install --no-deps -r "$MACE_PIP_LOCK"
else
  "$MM" run -n gmx-md python -m pip install -r "$ROOT/env/requirements-mace.txt"
fi
"$MM" run -n gmx-md python -m pip check
"$MM" run -n gmx-md python -c 'from mace.calculators import MACECalculator; import torch; print("MACE import OK; CPU/MPS/CUDA", True, torch.backends.mps.is_available(), torch.cuda.is_available())'
bash "$ROOT/scripts/setup/export_environment.sh"
