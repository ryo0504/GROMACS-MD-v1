#!/usr/bin/env bash
# CPU validation build, user-local prefix, independent of the gmx-md environment.
set -euo pipefail
source "$(dirname "$0")/common.sh"
[ "$(uname -s)" = Darwin ] || { echo 'Mac only; probe Linux first'; exit 2; }
for exe in "${GMX_BIN:-gmx}" gmx gmx_mpi; do
  if command -v "$exe" >/dev/null 2>&1; then
    echo "Existing GROMACS preserved: $(command -v "$exe")"; "$exe" --version; exit 0
  fi
done
if command -v brew >/dev/null 2>&1; then
  echo 'Available Homebrew route: brew install gromacs'
  echo 'To use this user-local source build instead: set GROMACS_LOCAL_BUILD=1'
  [ "${GROMACS_LOCAL_BUILD:-0}" = 1 ] || exit 0
fi
version="${GROMACS_VERSION:-2025.4}"
case "$version" in *[!0-9.]*|'') echo 'Invalid version'; exit 2;; esac
mkdir -p "$ROOT/.local/bootstrap"
if [ ! -x "$ROOT/.local/build-tools/bin/cmake" ]; then
  python3 -m venv "$ROOT/.local/build-tools"
  "$ROOT/.local/build-tools/bin/python" -m pip install 'cmake>=3.28,<4'
fi
curl -fL --retry 3 "https://ftp.gromacs.org/gromacs/gromacs-$version.tar.gz" -o "$ROOT/.local/bootstrap/gromacs-$version.tar.gz"
shasum -a 256 "$ROOT/.local/bootstrap/gromacs-$version.tar.gz" > "$ROOT/.local/bootstrap/gromacs-source.sha256"
# FFTW's Autotools/GNU make cannot reliably resolve macOS decomposed Unicode VPATH.
# An ASCII temporary build path also handles repositories with spaces.
build_root="$(mktemp -d /tmp/gmx-md-build.XXXXXX)"
printf '%s\n' "$build_root" > "$ROOT/.local/bootstrap/gromacs-build-directory.txt"
tar -xzf "$ROOT/.local/bootstrap/gromacs-$version.tar.gz" -C "$build_root"
cmake="$ROOT/.local/build-tools/bin/cmake"
"$cmake" -S "$build_root/gromacs-$version" -B "$build_root/build" \
  -DCMAKE_INSTALL_PREFIX="$ROOT/.local/gromacs" -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_GPU=OFF -DGMX_MPI=OFF -DGMX_OPENMP=OFF -DBUILD_TESTING=OFF -DGMX_BUILD_UNITTESTS=OFF -DCMAKE_BUILD_TYPE=Release
"$cmake" --build "$build_root/build" --parallel "${BUILD_JOBS:-4}"
"$cmake" --install "$build_root/build"
"$ROOT/.local/gromacs/bin/gmx" --version
