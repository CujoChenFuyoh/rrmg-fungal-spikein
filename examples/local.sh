#!/usr/bin/env bash
# Activate fungamr_meta first; this script never initializes Conda or calls PBS.
# Usage: RRMG_CONFIG=/absolute/path/config.sh bash examples/local.sh step2 1-100
set -euo pipefail
CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec bash "$CODE/rrmg_submit.sh" "${1:?Specify stage}" "${2:?Specify row or range}" --jobs "${LOCAL_JOBS:-1}"
