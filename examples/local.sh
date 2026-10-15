#!/usr/bin/env bash
# Activate fungamr_meta first; this script never initializes Conda or calls PBS.
# Usage: RRMG_CONFIG=/absolute/path/config.sh bash examples/local.sh step2 1-100 -t 1 -m 4
set -euo pipefail
CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAGE="${1:?Specify stage}"
ROWS="${2:?Specify row or range}"
shift 2
exec bash "$CODE/rrmg_submit.sh" "$STAGE" "$ROWS" --jobs "${LOCAL_JOBS:-1}" "$@"
