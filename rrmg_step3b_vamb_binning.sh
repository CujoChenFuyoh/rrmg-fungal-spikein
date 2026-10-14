#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  rrmg_step3b_vamb_binning.sh PLAN.tsv ROW_INDEX -w /path/to/project [options]

Purpose:
  Run one independent, single-sample VAMB binning analysis for one permutation:
    1. locate the run's QC reads and metaSPAdes contigs;
    2. map reads back to contigs and create a query-name-sorted BAM;
    3. run VAMB using contig composition and the single-sample abundance;
    4. align all contigs to the known spike-in genome;
    5. report per-bin and per-run target recovery metrics.

  This intentionally does not combine seeds in VAMB multi-sample mode because
  the current split RRMG design reuses the same background reads across seeds.

Required:
  PLAN.tsv and ROW_INDEX   Same permutation plan/index used by Steps 2 and 3
  -w DIR                   RRMG project working directory

Options:
  --target-ref FILE        Override target_fasta_resolved/target_fasta in plan
  -t, --threads INT        Threads for minimap2/samtools/VAMB; default 32
  --min-contig INT         Minimum contig length passed to VAMB; default 2000
  --min-bin-size INT       Minimum FASTA bin size emitted by VAMB; default 200000
  --vamb-seed INT          Fixed VAMB model seed; default 42
  --vamb-command FILE      VAMB executable name/path; default vamb
  --cuda                   Pass --cuda to VAMB
  --force                  Replace an incomplete/existing output for this run
  -h, --help               Show this help

Environment:
  Activate the analysis environment first; tools are resolved from PATH.

Outputs under each RUN_DIR:
  binning/read_to_contigs/*.name_sorted.bam
  binning/vamb/                         VAMB native output and FASTA bins
  binning/contigs.to_target.paf
  binning/vamb_bin_metrics.tsv
  binning/vamb_summary.tsv
  binning/VAMB_SUCCESS

Compatibility:
  VAMB 5.x uses "vamb bin default ... --bamdir".
  VAMB 3/4 uses "vamb ... --bamfiles".
  This wrapper detects and supports both command-line interfaces.
USAGE
}

PLAN="${1:-}"
ROW_INDEX="${2:-${PBS_ARRAY_INDEX:-}}"
if [[ -z "$PLAN" || -z "$ROW_INDEX" ]]; then
  usage
  exit 1
fi
shift 2 || true

WD=""
TARGET_REF_OVERRIDE=""
THREADS=32
MIN_CONTIG=2000
MIN_BIN_SIZE=200000
VAMB_SEED=42
VAMB_COMMAND="vamb"
USE_CUDA=0
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w) WD="$2"; shift 2;;
    --target-ref) TARGET_REF_OVERRIDE="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    --min-contig) MIN_CONTIG="$2"; shift 2;;
    --min-bin-size) MIN_BIN_SIZE="$2"; shift 2;;
    --vamb-seed) VAMB_SEED="$2"; shift 2;;
    --vamb-command) VAMB_COMMAND="$2"; shift 2;;
    --cuda) USE_CUDA=1; shift 1;;
    --force) FORCE=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 1;;
  esac
done

: "${WD:?Missing -w}"
for value_name in ROW_INDEX THREADS MIN_CONTIG MIN_BIN_SIZE VAMB_SEED; do
  value="${!value_name}"
  [[ "$value" =~ ^[0-9]+$ ]] || {
    echo "[ERROR] ${value_name} must be a non-negative integer: $value" >&2
    exit 1
  }
done
[[ "$ROW_INDEX" -ge 1 ]] || { echo "[ERROR] ROW_INDEX must be >= 1" >&2; exit 1; }
[[ "$THREADS" -ge 1 ]] || { echo "[ERROR] THREADS must be >= 1" >&2; exit 1; }
[[ "$MIN_CONTIG" -ge 1 ]] || { echo "[ERROR] MIN_CONTIG must be >= 1" >&2; exit 1; }
[[ "$MIN_BIN_SIZE" -ge 1 ]] || { echo "[ERROR] MIN_BIN_SIZE must be >= 1" >&2; exit 1; }

WD="$(mkdir -p "$WD" && cd "$WD" && pwd)"
cd "$WD"
[[ -s "$PLAN" ]] || { echo "[ERROR] Plan not found: $PLAN" >&2; exit 1; }
PLAN="$(cd "$(dirname "$PLAN")" && pwd)/$(basename "$PLAN")"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVALUATOR="${SCRIPT_DIR}/rrmg_evaluate_vamb_bins.py"
[[ -s "$EVALUATOR" ]] || {
  echo "[ERROR] Evaluator not found: $EVALUATOR" >&2
  exit 1
}

# All programs use the caller's active PATH, including binary BAM streams.

mkdir -p "${WD}/temp"
ROW_ENV="$(mktemp "${WD}/temp/rrmg_vamb_row.XXXXXX")"
python3 - "$PLAN" "$ROW_INDEX" "$WD" > "$ROW_ENV" <<'PY'
import csv
import re
import shlex
import sys
from pathlib import Path

plan = Path(sys.argv[1])
row_index = int(sys.argv[2])
wd = Path(sys.argv[3])

def safe(value):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return value.strip("_") or "NA"

def q(value):
    return shlex.quote(str(value))

chosen = None
with plan.open(newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if int(row.get("run_index", "0") or 0) == row_index:
            chosen = row
            break
if chosen is None:
    raise SystemExit("[ERROR] Plan row not found: {}".format(row_index))

sample = chosen.get("sample") or chosen.get("run_id") or "row{}".format(row_index)
seed = chosen.get("seed") or row_index
run_id = chosen.get("run_id") or "{}__Seed{}".format(sample, seed)
run_dir = wd / "temp" / "rrmg_permutation" / safe(run_id)
target_ref = (
    chosen.get("target_fasta_resolved")
    or chosen.get("target_fasta")
    or ""
)
values = {
    "RUN_INDEX": row_index,
    "RUN_ID": run_id,
    "SAMPLE": sample,
    "SEED": seed,
    "RUN_DIR": run_dir,
    "QC_R1": run_dir / "qc" / "{}.Seed{}_1.fastq.gz".format(sample, seed),
    "QC_R2": run_dir / "qc" / "{}.Seed{}_2.fastq.gz".format(sample, seed),
    "MIX_R1": run_dir / "reads" / "{}.Seed{}_1.fastq.gz".format(safe(sample), seed),
    "MIX_R2": run_dir / "reads" / "{}.Seed{}_2.fastq.gz".format(safe(sample), seed),
    "CONTIGS": run_dir / "asm" / "contigs.fasta",
    "TARGET_REF_FROM_PLAN": target_ref,
    "TARGET_ABUNDANCE": chosen.get("target_abundance", ""),
    "TOTAL_PAIRS": chosen.get("total_pairs", ""),
    "TARGET_EXPECTED_COV": chosen.get("target_expected_cov", ""),
}
for key, value in values.items():
    print("{}={}".format(key, q(value)))
PY
# shellcheck disable=SC1090
source "$ROW_ENV"
rm -f "$ROW_ENV"

TARGET_REF="${TARGET_REF_OVERRIDE:-$TARGET_REF_FROM_PLAN}"
if [[ -n "$TARGET_REF" && "$TARGET_REF" != /* ]]; then
  TARGET_REF="${WD}/${TARGET_REF}"
fi

[[ -s "$CONTIGS" ]] || {
  echo "[ERROR] Assembly contigs not found; complete Step 3 first: $CONTIGS" >&2
  exit 1
}
[[ -s "$TARGET_REF" ]] || {
  echo "[ERROR] Target reference not found. Set --target-ref or target_fasta_resolved in plan: $TARGET_REF" >&2
  exit 1
}

READ_R1="$QC_R1"
READ_R2="$QC_R2"
READ_SOURCE="qc"
if [[ ! -s "$READ_R1" || ! -s "$READ_R2" ]]; then
  READ_R1="$MIX_R1"
  READ_R2="$MIX_R2"
  READ_SOURCE="mixed"
fi
[[ -s "$READ_R1" && -s "$READ_R2" ]] || {
  echo "[ERROR] Neither QC nor mixed read pair was found for $RUN_ID" >&2
  exit 1
}

BINNING_DIR="${RUN_DIR}/binning"
BAM_DIR="${BINNING_DIR}/read_to_contigs"
VAMB_OUT="${BINNING_DIR}/vamb"
BAM="${BAM_DIR}/${SAMPLE}.Seed${SEED}.to_contigs.name_sorted.bam"
MMI="${BAM_DIR}/contigs.mmi"
FLAGSTAT="${BAM}.flagstat.txt"
MAP_LOG="${BAM}.minimap2.log"
TARGET_PAF="${BINNING_DIR}/contigs.to_target.paf"
TARGET_MAP_LOG="${TARGET_PAF}.minimap2.log"
PER_BIN_TSV="${BINNING_DIR}/vamb_bin_metrics.tsv"
SUMMARY_TSV="${BINNING_DIR}/vamb_summary.tsv"
SUCCESS_MARKER="${BINNING_DIR}/VAMB_SUCCESS"

mkdir -p "$BINNING_DIR" "$BAM_DIR"

if [[ "$FORCE" -eq 0 && -s "$SUCCESS_MARKER" && -s "$SUMMARY_TSV" ]]; then
  echo "[SKIP] Completed VAMB result already exists: $SUMMARY_TSV"
  exit 0
fi

if [[ -e "$VAMB_OUT" ]]; then
  if [[ "$FORCE" -eq 1 ]]; then
    case "$VAMB_OUT" in
      "$RUN_DIR"/binning/vamb) rm -rf -- "$VAMB_OUT";;
      *) echo "[ERROR] Refusing to remove unexpected VAMB path: $VAMB_OUT" >&2; exit 1;;
    esac
  else
    echo "[ERROR] Existing incomplete VAMB output: $VAMB_OUT" >&2
    echo "[ERROR] Inspect it, then rerun this one row with --force if replacement is intended." >&2
    exit 1
  fi
fi
rm -f "$SUCCESS_MARKER"

for executable in minimap2 samtools python3; do
  if ! "$executable" --version >/dev/null 2>&1; then
    if [[ "$executable" == "python3" ]] && python3 -V >/dev/null 2>&1; then
      continue
    fi
    echo "[ERROR] Command is unavailable in conda env ${RRMG_CONDA_ENV}: $executable" >&2
    exit 1
  fi
done
if ! "$VAMB_COMMAND" --help >/dev/null 2>&1; then
  echo "[ERROR] VAMB is unavailable in conda env ${RRMG_CONDA_ENV}: $VAMB_COMMAND" >&2
  exit 1
fi

VAMB_VERSION="$("$VAMB_COMMAND" --version 2>&1 | head -n 1 || true)"
VAMB_MODE=""
if "$VAMB_COMMAND" bin default --help >/dev/null 2>&1; then
  VAMB_MODE="v5"
elif "$VAMB_COMMAND" --help 2>&1 | grep -q -- "--bamfiles"; then
  VAMB_MODE="legacy"
else
  echo "[ERROR] Unsupported VAMB command-line interface: ${VAMB_COMMAND}" >&2
  exit 1
fi

echo "++++++++++++++++++++++++++++++++++++++++"
echo "RRMG Step 3B: one-permutation VAMB binning"
echo "Run index       : $RUN_INDEX"
echo "Run ID          : $RUN_ID"
echo "Reads           : $READ_SOURCE ($READ_R1 / $READ_R2)"
echo "Contigs         : $CONTIGS"
echo "Target reference: $TARGET_REF"
echo "Conda env       : $RRMG_CONDA_ENV"
echo "VAMB mode       : $VAMB_MODE"
echo "VAMB version    : $VAMB_VERSION"
echo "VAMB model seed : $VAMB_SEED"
echo "Threads         : $THREADS"
echo "++++++++++++++++++++++++++++++++++++++++"

minimap2 -I 32G -d "$MMI" "$CONTIGS" 2> "$MAP_LOG"

# VAMB 3.x requires alignments for one read to remain consecutive.  A
# query-name-sorted BAM satisfies that requirement and also works with VAMB 5.
minimap2 -I 32G -t "$THREADS" -ax sr "$MMI" "$READ_R1" "$READ_R2" \
  2>> "$MAP_LOG" \
  | samtools view -@ "$THREADS" -F 3584 -b - \
  | samtools sort -n -@ "$THREADS" -m 1G -O BAM -o "$BAM" -

samtools quickcheck -v "$BAM"
samtools flagstat -@ "$THREADS" "$BAM" > "$FLAGSTAT"

VAMB_ARGS=(
  --outdir "$VAMB_OUT"
  --fasta "$CONTIGS"
  -p "$THREADS"
  -m "$MIN_CONTIG"
  --minfasta "$MIN_BIN_SIZE"
  --seed "$VAMB_SEED"
)
if [[ "$USE_CUDA" -eq 1 ]]; then
  VAMB_ARGS+=(--cuda)
fi

if [[ "$VAMB_MODE" == "v5" ]]; then
  "$VAMB_COMMAND" bin default "${VAMB_ARGS[@]}" --bamdir "$BAM_DIR"
else
  "$VAMB_COMMAND" "${VAMB_ARGS[@]}" --bamfiles "$BAM"
fi

CLUSTERS_TSV=""
for candidate in \
  "${VAMB_OUT}/vae_clusters_unsplit.tsv" \
  "${VAMB_OUT}/vae_clusters_split.tsv" \
  "${VAMB_OUT}/clusters.tsv"; do
  if [[ -s "$candidate" ]]; then
    CLUSTERS_TSV="$candidate"
    break
  fi
done
[[ -n "$CLUSTERS_TSV" ]] || {
  echo "[ERROR] VAMB completed without a recognized cluster TSV in $VAMB_OUT" >&2
  exit 1
}

BINS_DIR="${VAMB_OUT}/bins"
minimap2 -t "$THREADS" -x asm5 "$TARGET_REF" "$CONTIGS" \
  > "$TARGET_PAF" 2> "$TARGET_MAP_LOG"

python3 "$EVALUATOR" \
  --clusters "$CLUSTERS_TSV" \
  --contigs "$CONTIGS" \
  --target-ref "$TARGET_REF" \
  --target-paf "$TARGET_PAF" \
  --bins-dir "$BINS_DIR" \
  --out-bins "$PER_BIN_TSV" \
  --out-summary "$SUMMARY_TSV" \
  --run-index "$RUN_INDEX" \
  --run-id "$RUN_ID" \
  --sample "$SAMPLE" \
  --seed "$SEED" \
  --target-abundance "$TARGET_ABUNDANCE" \
  --total-pairs "$TOTAL_PAIRS" \
  --target-expected-cov "$TARGET_EXPECTED_COV" \
  --vamb-mode "$VAMB_MODE" \
  --vamb-version "$VAMB_VERSION" \
  --vamb-seed "$VAMB_SEED" \
  --min-contig "$MIN_CONTIG" \
  --min-bin-size "$MIN_BIN_SIZE" \
  --bam "$BAM" \
  --vamb-outdir "$VAMB_OUT"

printf 'run_id=%s\ncompleted=%s\n' "$RUN_ID" "$(date -Is)" > "$SUCCESS_MARKER"
echo "[DONE] VAMB summary: $SUMMARY_TSV"
