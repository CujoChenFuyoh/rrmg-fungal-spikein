#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  rrmg_submit.sh STEP [ROW|START-END] [options]
  rrmg_run_local.sh STEP [ROW|START-END] [options]  # compatible alias

STEP:
  step1   background-read simulation
  step2   target-read simulation and merge
  step3   QC, optional Kraken2, assembly, and target-gene mapping
  step3b  one-permutation VAMB binning
  step4   QC-read to target-gene depth
  step5   target-only read to target-genome coverage
  step45  optional read-depth and target-genome coverage analyses

Options:
  -t, --threads N   Threads per row; overrides the selected step's config value
  -m, --memory GB   Per-row memory in GiB for step1/2 Java or step3 SPAdes
                   (positive integer, no workflow-imposed maximum)
  -j, --jobs N      Run at most N rows concurrently (default: LOCAL_JOBS or 1)
  --force           Override the selected step's force setting for this run
  --dry-run         Validate inputs and print rows without running them
  -h, --help        Show this help

Examples:
  ./rrmg_run_local.sh step1
  ./rrmg_run_local.sh step2 215-480
  ./rrmg_run_local.sh step2 215-480 -t 1 -m 4 --jobs 4
  ./rrmg_run_local.sh step3 215 --threads 16 --memory 64 --force

Activate the analysis environment first. This runner never uses qsub or .pbs
files. ROW is run_index (bg_index for step1), not Seed or a text line number.
--memory controls supported application settings; PBS/Slurm memory must still be
requested in the scheduler script. step3b/4/5/45 have no total-memory CLI cap.
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
CLI_THREADS=""
CLI_MEMORY_GB=""
CLI_JOBS=""
FORCE=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--threads)
      [[ $# -ge 2 ]] || { echo "[ERROR] $1 requires a value" >&2; exit 1; }
      CLI_THREADS="$2"
      shift 2
      ;;
    -m|--memory)
      [[ $# -ge 2 ]] || { echo "[ERROR] $1 requires a value" >&2; exit 1; }
      CLI_MEMORY_GB="$2"
      shift 2
      ;;
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
if [[ -n "$CLI_THREADS" && ! "$CLI_THREADS" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] --threads must be a positive integer: $CLI_THREADS" >&2
  exit 1
fi
if [[ -n "$CLI_MEMORY_GB" && ! "$CLI_MEMORY_GB" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] --memory must be a positive integer in GiB: $CLI_MEMORY_GB" >&2
  exit 1
fi
if [[ -n "$CLI_MEMORY_GB" ]]; then
  case "$STEP" in
    step1|step2|step3) ;;
    *)
      echo "[ERROR] --memory has no application-level meaning for $STEP." >&2
      echo "        Request memory in the scheduler script; omit --memory for this step." >&2
      exit 1
      ;;
  esac
fi
if [[ "$STEP" == step2 && -n "$CLI_THREADS" && "$CLI_THREADS" -ne 1 ]]; then
  echo "[WARN] Step2 uses one thread by default for fixed-seed reproducibility." >&2
fi

case "$STEP" in
  step1)
    TABLE="${BACKGROUND_PLAN:-temp/rrmg_background/background_plan.tsv}"
    ;;
  step2)
    TABLE="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
    ;;
  step3)
    TABLE="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
    ;;
  step3b)
    TABLE="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
    ;;
  step4|step5|step45)
    TABLE="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
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

WRAPPER_PATH="$PIPELINE_DIR/rrmg_run_one.sh"
[[ -s "$TABLE" ]] || { echo "[ERROR] Input table not found: $TABLE" >&2; exit 1; }
[[ -s "$WRAPPER_PATH" ]] || { echo "[ERROR] Stage wrapper not found: $WRAPPER_PATH" >&2; exit 1; }

# Validate explicit IDs before launching anything. Never renumber a subset.
ROWS_TEXT="$(python3 - "$TABLE" "$STEP" "$RANGE" <<'PY'
import csv, re, sys
table, step, selection = sys.argv[1:]
key = 'bg_index' if step == 'step1' else 'run_index'
with open(table, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    if key not in (reader.fieldnames or []):
        raise SystemExit('[ERROR] Missing index column: ' + key)
    raw = [r[key] for r in reader]
if not raw or any(not re.fullmatch(r'[1-9][0-9]*', s or '') for s in raw):
    raise SystemExit('[ERROR] Plan indices must be positive integers')
ids = list(map(int, raw))
if len(ids) != len(set(ids)):
    raise SystemExit('[ERROR] Duplicate plan indices')
chosen = ids
if selection:
    if not re.fullmatch(r'[1-9][0-9]*(?:-[1-9][0-9]*)?', selection):
        raise SystemExit('[ERROR] Use ROW or START-END')
    bounds = list(map(int, selection.split('-')))
    start, end = bounds[0], bounds[-1]
    if end < start or end - start + 1 > len(ids):
        raise SystemExit('[ERROR] Invalid row selection: ' + selection)
    chosen = list(range(start, end + 1))
    missing = sorted(set(chosen) - set(ids))
    if missing:
        raise SystemExit('[ERROR] Indices absent from plan: ' + str(missing))
# Write explicit LF bytes so Git Bash never receives a stray CR when a Windows
# Python executable is used for local validation.
sys.stdout.buffer.write(('\n'.join(map(str, chosen))).encode())
PY
)"
mapfile -t ROWS <<< "$ROWS_TEXT"

RUN_ID="${RRMG_LOCAL_RUN_ID:-local_${STEP}_$(date +%Y%m%d_%H%M%S)_$$}"
LOG_DIR="${RRMG_LOCAL_LOG_DIR:-$WD/logs/local/$RUN_ID}"
mkdir -p "$LOG_DIR"
SUMMARY="$LOG_DIR/summary.tsv"
printf 'row_index\tstatus\n' > "$SUMMARY"

echo "[LOCAL] step=$STEP selection=${RANGE:-all} rows=${#ROWS[@]} jobs=$JOBS threads=${CLI_THREADS:-from-config} memory_gb=${CLI_MEMORY_GB:-from-config/not-applicable} force=${FORCE:-from-config}"
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
  RRMG_LOG_DIR="$LOG_DIR" \
  RRMG_THREADS_OVERRIDE="$CLI_THREADS" \
  RRMG_MEMORY_GB_OVERRIDE="$CLI_MEMORY_GB" \
  RRMG_FORCE_OVERRIDE="$FORCE" \
    bash "$WRAPPER_PATH" "$STEP" "$row" > "$LOG_DIR/row_${row}.log" 2>&1
}

if (( DRY_RUN == 1 )); then
  for row in "${ROWS[@]}"; do
    printf '[DRY-RUN] RRMG_CONFIG=%q bash %q %q %q\n' "$CONFIG" "$WRAPPER_PATH" "$STEP" "$row"
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

for row in "${ROWS[@]}"; do
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
