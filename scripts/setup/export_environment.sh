#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
platform="$($MM info --json | python3 -c 'import json,sys; print(json.load(sys.stdin)["platform"])')"
"$MM" env export -n gmx-md --explicit > "$ROOT/env/locks/$platform.lock"
"$MM" run -n gmx-md python "$ROOT/scripts/setup/export_packages.py" "$ROOT/env/locks" "$platform"
