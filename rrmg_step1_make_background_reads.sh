#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  rrmg_step1_make_background_reads.sh \
    --background-plan temp/rrmg_background/background_plan.tsv \
    --row-index ${PBS_ARRAY_INDEX} \
    -w /path/to/project

Purpose:
  Generate ONE fixed background-read set from one row of background_plan.tsv.
  This is Step 1B in the batched background workflow.

Input:
  background_plan.tsv -> selected background_design.tsv -> background FASTAs

Output:
  <bg_dir>/reads/<bg_key>.background_1.fastq.gz
  <bg_dir>/reads/<bg_key>.background_2.fastq.gz

Run after:
  rrmg_make_background_plan.py

Required environment:
  RANDOMREADSMG_SH      Absolute path to randomreadsmg.sh

Optional environment:
  THREADS               RandomReadsMG t= value. Default: 8
  BBMAP_XMX             Java heap for BBMap wrapper. Default: 32g
  BACKGROUND_SEED       Fixed background seed. Default: 1000

Options:
  --force               Regenerate even if merged background reads already exist
USAGE
}

BACKGROUND_PLAN=""
ROW_INDEX="${PBS_ARRAY_INDEX:-}"
WD=""
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --background-plan) BACKGROUND_PLAN="$2"; shift 2;;
    --row-index) ROW_INDEX="$2"; shift 2;;
    -w) WD="$2"; shift 2;;
    --force) FORCE=1; shift 1;;
    --keep-per-genome) shift 1;; # retained for backward compatibility; no per-genome FASTQs are generated now
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

: "${BACKGROUND_PLAN:?Missing --background-plan}"
: "${ROW_INDEX:?Missing --row-index or PBS_ARRAY_INDEX}"
: "${WD:?Missing -w}"
: "${RANDOMREADSMG_SH:?RANDOMREADSMG_SH is not set; export absolute path to randomreadsmg.sh}"

THREADS="${THREADS:-8}"
BBMAP_XMX="${BBMAP_XMX:-32g}"
BACKGROUND_SEED="${BACKGROUND_SEED:-1000}"

WD="$(mkdir -p "$WD" && cd "$WD" && pwd)"
cd "$WD"

[[ -s "$BACKGROUND_PLAN" ]] || { echo "[ERROR] Background plan not found: $BACKGROUND_PLAN" >&2; exit 1; }
[[ -s "$RANDOMREADSMG_SH" ]] || { echo "[ERROR] RANDOMREADSMG_SH not found: $RANDOMREADSMG_SH" >&2; exit 1; }

mkdir -p "${WD}/temp"
ROW_ENV="$(mktemp "${WD}/temp/rrmg_bg_row_env.XXXXXX")"
python3 - "$BACKGROUND_PLAN" "$ROW_INDEX" > "$ROW_ENV" <<'PY'
import csv
import shlex
import sys
from pathlib import Path

plan = Path(sys.argv[1])
row_index = int(sys.argv[2])

with plan.open(newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    chosen = None
    for row in reader:
        idx = int(float(row.get("bg_index", "0") or 0))
        if idx == row_index:
            chosen = row
            break

if chosen is None:
    raise SystemExit("[ERROR] Background row not found: {}".format(row_index))

def q(v):
    return shlex.quote(str(v))

for k, v in chosen.items():
    print("{}={}".format(k.upper(), q(v)))
PY

# shellcheck disable=SC1090
source "$ROW_ENV"
rm -f "$ROW_ENV"

mkdir -p "${BG_DIR}/reads" "${BG_DIR}/logs"

log="${BG_DIR}/reads/${BG_KEY}.randomreadsmg.log"

cat <<INFO
[BG] bg_index=${BG_INDEX}
[BG] bg_key=${BG_KEY}
[BG] target_abundance=${TARGET_ABUNDANCE}
[BG] total_pairs=${TOTAL_PAIRS}
[BG] background_pairs=${BACKGROUND_PAIRS}
[BG] readlen=${READLEN}
[BG] design=${DESIGN_TSV}
[BG] bg_r1=${BG_R1}
[BG] bg_r2=${BG_R2}
[BG] THREADS=${THREADS}
[BG] BBMAP_XMX=${BBMAP_XMX}
[BG] BACKGROUND_SEED=${BACKGROUND_SEED}
INFO

[[ -s "$DESIGN_TSV" ]] || { echo "[ERROR] Missing design TSV: $DESIGN_TSV" >&2; exit 1; }

if [[ "$FORCE" -eq 0 && -s "$BG_R1" && -s "$BG_R2" ]]; then
  echo "[BG] background exists, skip: $BG_R1 / $BG_R2"
  exit 0
fi

MG_ARGS=()
while IFS=$'\t' read -r idx gid name taxid role rel_final bg_weight fasta glen sim_pairs sim_cov sim_seed; do
  [[ "$idx" == "idx" ]] && continue
  [[ -s "$fasta" ]] || { echo "[ERROR] FASTA not found for ${gid}: $fasta" >&2; exit 1; }

  # RandomReadsMG custom coverage syntax: /path/genome.fna=0.123
  # sim_cov is already calculated in background_design.tsv.
  MG_ARGS+=("${fasta}=${sim_cov}")
  echo "[BG-ADD] ${gid} pairs=${sim_pairs} cov=${sim_cov} fasta=${fasta}"
done < "$DESIGN_TSV"

[[ "${#MG_ARGS[@]}" -gt 0 ]] || {
  echo "[ERROR] No background genomes found in design: $DESIGN_TSV" >&2
  exit 1
}

rm -f "$BG_R1" "$BG_R2"

bash "$RANDOMREADSMG_SH" \
  -Xmx"$BBMAP_XMX" \
  "${MG_ARGS[@]}" \
  out="$BG_R1" \
  out2="$BG_R2" \
  len="$READLEN" \
  paired=t \
  t="$THREADS" \
  seed="$BACKGROUND_SEED" \
  zl=4 \
  2> "$log"

[[ -s "$BG_R1" && -s "$BG_R2" ]] || {
  echo "[ERROR] background randomreadsmg failed; see $log" >&2
  exit 1
}

echo "[BG-DONE] ${BG_KEY}"
echo "[BG-DONE] ${BG_R1}"
echo "[BG-DONE] ${BG_R2}"
