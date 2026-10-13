#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  rrmg_run_local.sh STEP [ROW|START-END] [options]

STEP:
  step1   background-read simulation
  step2   target-read simulation and merge
  step3   QC, optional Kraken2, assembly, and target-gene mapping
  step3b  one-permutation VAMB binning
  step45  optional read-depth and target-genome coverage analyses

Options:
  -j, --jobs N   Run at most N rows concurrently (default: LOCAL_JOBS or 1)
  --force        Override the selected step's force setting for this run
  --dry-run      Validate inputs and print rows without running them
  -h, --help     Show this help

Examples:
  ./rrmg_run_local.sh step1
  ./rrmg_run_local.sh step2 215-480
  ./rrmg_run_local.sh step2 215-480 --jobs 4
  ./rrmg_run_local.sh step3 215 --force

This runner does not use qsub. It executes the maintained stage wrappers
directly, so PBS and non-PBS runs share the same configuration and arguments.
USAGE
}

[[ $# -ge 1 ]] || { usage; exit 1; }

STEP="$1"
shift
if [[ "$STEP" == "-h" || "$STEP" == "--help" ]]; then
  usage
  exit 0
fi

RANGE=""
CLI_JOBS=""
FORCE=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -j|--jobs)
      [[ $# -ge 2 ]] || { echo "[ERROR] $1 requires a value" >&2; exit 1; }
      CLI_JOBS="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    [0-9]*-[0-9]*|[0-9]*)
      [[ -z "$RANGE" ]] || { echo "[ERROR] More than one row range was supplied" >&2; exit 1; }
      RANGE="$1"
      shift
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${RRMG_CONFIG:-${SCRIPT_DIR}/rrmg_config.sh}"
[[ -r "$CONFIG" ]] || {
  echo "[ERROR] Configuration not found: $CONFIG" >&2
  echo "        cp ${SCRIPT_DIR}/rrmg_config.example.sh ${SCRIPT_DIR}/rrmg_config.sh" >&2
  exit 1
}
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"

# shellcheck disable=SC1090
source "$CONFIG"
: "${PIPELINE_DIR:?Set PIPELINE_DIR in $CONFIG}"
: "${WD:?Set WD in $CONFIG}"

JOBS="${CLI_JOBS:-${LOCAL_JOBS:-1}}"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || {
  echo "[ERROR] --jobs must be a positive integer: $JOBS" >&2
  exit 1
}

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

WRAPPER_PATH="$PIPELINE_DIR/$WRAPPER"
[[ -s "$TABLE" ]] || { echo "[ERROR] Input table not found: $TABLE" >&2; exit 1; }
[[ -s "$WRAPPER_PATH" ]] || { echo "[ERROR] Stage wrapper not found: $WRAPPER_PATH" >&2; exit 1; }

ROW_COUNT="$(awk 'NR > 1 {n++} END {print n+0}' "$TABLE")"
[[ "$ROW_COUNT" -ge 1 ]] || { echo "[ERROR] No data rows in: $TABLE" >&2; exit 1; }

if [[ -z "$RANGE" ]]; then
  START=1
  END="$ROW_COUNT"
elif [[ "$RANGE" =~ ^([0-9]+)-([0-9]+)$ ]]; then
  START="${BASH_REMATCH[1]}"
  END="${BASH_REMATCH[2]}"
elif [[ "$RANGE" =~ ^[0-9]+$ ]]; then
  START="$RANGE"
  END="$RANGE"
else
  echo "[ERROR] Invalid row selection: $RANGE" >&2
  exit 1
fi

if (( START < 1 || END < START || END > ROW_COUNT )); then
  echo "[ERROR] Invalid row selection ${START}-${END}; valid rows are 1-${ROW_COUNT}" >&2
  exit 1
fi

RUN_ID="${RRMG_LOCAL_RUN_ID:-local_${STEP}_$(date +%Y%m%d_%H%M%S)_$$}"
LOG_DIR="${RRMG_LOCAL_LOG_DIR:-$WD/logs/local/$RUN_ID}"
mkdir -p "$LOG_DIR"
SUMMARY="$LOG_DIR/summary.tsv"
printf 'row_index\tstatus\n' > "$SUMMARY"

echo "[LOCAL] step=$STEP rows=${START}-${END} jobs=$JOBS force=$FORCE"
echo "[LOCAL] table=$TABLE"
echo "[LOCAL] logs=$LOG_DIR"
if (( JOBS > 1 )); then
  echo "[WARN] Concurrent rows multiply the per-row CPU and memory requirements."
fi

run_row() {
  local row="$1"
  echo "[START] step=$STEP row=$row"
  RRMG_CONFIG="$CONFIG" \
  RRMG_ROW_INDEX="$row" \
  PBS_ARRAY_INDEX="$row" \
  PBS_JOBID="$RUN_ID" \
  RRMG_LOG_DIR="$LOG_DIR" \
  RRMG_FORCE_OVERRIDE="$FORCE" \
    bash "$WRAPPER_PATH"
}

if (( DRY_RUN == 1 )); then
  for ((row=START; row<=END; row++)); do
    printf '[DRY-RUN] row=%d RRMG_CONFIG=%q bash %q\n' "$row" "$CONFIG" "$WRAPPER_PATH"
    printf '%s\tDRY_RUN\n' "$row" >> "$SUMMARY"
  done
  echo "[DONE] Dry run completed: $SUMMARY"
  exit 0
fi

declare -a PIDS=()
declare -a BATCH_ROWS=()
declare -a FAILED_ROWS=()
SUCCEEDED=0

stop_children() {
  trap - INT TERM
  echo "[WARN] Interrupt received; stopping active local rows" >&2
  for pid in "${PIDS[@]:-}"; do
    kill "$pid" 2>/dev/null || true
  done
  wait 2>/dev/null || true
  exit 130
}
trap stop_children INT TERM

reap_batch() {
  local i pid row
  for i in "${!PIDS[@]}"; do
    pid="${PIDS[$i]}"
    row="${BATCH_ROWS[$i]}"
    if wait "$pid"; then
      echo "[OK] step=$STEP row=$row"
      printf '%s\tOK\n' "$row" >> "$SUMMARY"
      SUCCEEDED=$((SUCCEEDED + 1))
    else
      echo "[FAIL] step=$STEP row=$row" >&2
      printf '%s\tFAILED\n' "$row" >> "$SUMMARY"
      FAILED_ROWS+=("$row")
    fi
  done
  PIDS=()
  BATCH_ROWS=()
}

for ((row=START; row<=END; row++)); do
  run_row "$row" &
  PIDS+=("$!")
  BATCH_ROWS+=("$row")
  if (( ${#PIDS[@]} >= JOBS )); then
    reap_batch
  fi
done

if (( ${#PIDS[@]} > 0 )); then
  reap_batch
fi

FAILED_COUNT="${#FAILED_ROWS[@]}"
echo "[SUMMARY] succeeded=$SUCCEEDED failed=$FAILED_COUNT file=$SUMMARY"
if (( FAILED_COUNT > 0 )); then
  printf '[ERROR] Failed rows:' >&2
  printf ' %s' "${FAILED_ROWS[@]}" >&2
  printf '\n' >&2
  exit 1
fi

echo "[DONE] Local workflow rows completed successfully"
