#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  rrmg_submit.sh STEP [START-END] [extra qsub options]

STEP:
  step1   fixed background-read simulation
  step2   target-read simulation and merge
  step3   QC, optional Kraken2, assembly, and target-gene mapping
  step3b  one-permutation VAMB binning
  step45  optional read-depth and target-genome coverage analyses

Examples:
  ./rrmg_submit.sh step1
  ./rrmg_submit.sh step2 215-480
  ./rrmg_submit.sh step3 1-300 -W depend=afterok:12345

Configuration:
  Set RRMG_CONFIG or copy rrmg_config.example.sh to rrmg_config.sh.
USAGE
}

[[ $# -ge 1 ]] || { usage; exit 1; }

STEP="$1"
shift

if [[ "$STEP" == "-h" || "$STEP" == "--help" ]]; then
  usage
  exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${RRMG_CONFIG:-${SCRIPT_DIR}/rrmg_config.sh}"
[[ -r "$CONFIG" ]] || {
  echo "[ERROR] Configuration not found: $CONFIG" >&2
  exit 1
}
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"

# shellcheck disable=SC1090
source "$CONFIG"
: "${PIPELINE_DIR:?Set PIPELINE_DIR in $CONFIG}"
: "${WD:?Set WD in $CONFIG}"

case "$STEP" in
  step1)
    TABLE="${BACKGROUND_PLAN:-temp/rrmg_background/background_plan.tsv}"
    WRAPPER="rrmg_step1_background_array.pbs"
    ;;
  step2)
    TABLE="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
    WRAPPER="rrmg_step2_target_merge_array.pbs"
    ;;
  step3)
    TABLE="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
    WRAPPER="rrmg_step3_downstream_array.pbs"
    ;;
  step3b)
    TABLE="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
    WRAPPER="rrmg_step3b_vamb_array.pbs"
    ;;
  step45)
    TABLE="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
    WRAPPER="rrmg_step4_5_supplementary_downstream_array.pbs"
    ;;
  *)
    echo "[ERROR] Unknown step: $STEP" >&2
    usage
    exit 1
    ;;
esac

case "$TABLE" in
  /*) ;;
  *) TABLE="$WD/$TABLE" ;;
esac

[[ -s "$TABLE" ]] || { echo "[ERROR] Input table not found: $TABLE" >&2; exit 1; }
[[ -s "$PIPELINE_DIR/$WRAPPER" ]] || { echo "[ERROR] PBS wrapper not found: $PIPELINE_DIR/$WRAPPER" >&2; exit 1; }
command -v qsub >/dev/null 2>&1 || { echo "[ERROR] qsub is not available" >&2; exit 1; }

ROW_COUNT="$(awk 'NR > 1 {n++} END {print n+0}' "$TABLE")"
[[ "$ROW_COUNT" -ge 1 ]] || { echo "[ERROR] No data rows in: $TABLE" >&2; exit 1; }

RANGE="1-${ROW_COUNT}"
if [[ $# -ge 1 && "$1" =~ ^[0-9]+-[0-9]+$ ]]; then
  RANGE="$1"
  shift
fi

START="${RANGE%-*}"
END="${RANGE#*-}"
if (( START < 1 || END < START || END > ROW_COUNT )); then
  echo "[ERROR] Invalid array range ${RANGE}; valid rows are 1-${ROW_COUNT}" >&2
  exit 1
fi

echo "[SUBMIT] step=$STEP rows=$RANGE table=$TABLE"
qsub -J "$RANGE" -v "RRMG_CONFIG=$CONFIG" "$@" "$PIPELINE_DIR/$WRAPPER"
