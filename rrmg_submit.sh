#!/usr/bin/env bash
# Historical name, scheduler-independent since 0.3.0. Never calls qsub.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$SCRIPT_DIR/rrmg_run_local.sh" "$@"
