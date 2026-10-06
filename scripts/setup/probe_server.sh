#!/usr/bin/env bash
# Read-only, tolerates unavailable commands. Does not load modules or alter drivers.
set -uo pipefail
source "$(dirname "$0")/common.sh"
report="${1:-$ROOT/cluster/detected/platform_report.txt}"
mkdir -p "$(dirname "$report")"
probe() {
  printf '\n## %s\n' "$*"
  if command -v "$1" >/dev/null 2>&1; then "$@" 2>&1 || true; else echo '未導入'; fi
}
{
  date -u
  probe uname -a
  probe uname -m
  if [ -r /etc/os-release ]; then cat /etc/os-release; fi
  probe lscpu
  probe sysctl hw.physicalcpu hw.logicalcpu hw.memsize machdep.cpu.brand_string
  probe free -h
  probe nvidia-smi
  probe nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
  probe nvcc --version
  probe gcc --version
  probe clang --version
  probe cmake --version
  probe mpirun --version
  probe mpicc --version
  probe ompi_info --version
  probe sbatch --version
  probe srun --version
  probe sinfo
  if command -v module >/dev/null 2>&1; then
    probe module list
    probe module avail
  else
    echo 'module: unavailable in this non-login shell; also run module list/avail in your cluster login shell.'
  fi
  probe gmx --version
  probe gmx_mpi --version
  probe df -h "$ROOT" "$HOME"
  if [ -n "${SCRATCH:-}" ] && [ -d "$SCRATCH" ]; then probe df -h "$SCRATCH"; fi
  probe mount
  probe quota -s
  printf '\nCPU/GPU/MPI/SLURM choices require cluster policy and this report.\n'
} > "$report" 2>&1
echo "$report"
