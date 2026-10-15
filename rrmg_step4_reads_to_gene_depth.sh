#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTQ_MIXER="${FASTQ_MIXER:-$SCRIPT_DIR/rrmg_mix_paired_fastq.py}"

usage() {
  cat <<'USAGE'
Usage:
  rrmg_step4_reads_to_gene_depth.sh plan.tsv ROW_INDEX \
    -w /path/to/project \
    -gene ref/CYP51A_wt.fasta

Purpose:
  For one RRMG permutation run, map QC or raw paired reads to the target gene
  and calculate per-base read depth over the gene.

Inputs expected in each run directory:
  temp/rrmg_permutation/<run_id>/qc/*_1.fastq.gz
  temp/rrmg_permutation/<run_id>/qc/*_2.fastq.gz

Outputs:
  temp/rrmg_permutation/<run_id>/readmap/reads_to_GENE.bam
  temp/rrmg_permutation/<run_id>/readmap/reads_to_GENE.bam.bai
  temp/rrmg_permutation/<run_id>/readmap/reads_to_GENE.depth.tsv
  temp/rrmg_permutation/<run_id>/readmap/reads_to_GENE.coverage_summary.tsv
  temp/rrmg_permutation/<run_id>/readmap/reads_to_GENE.paf.gz

Required tools:
  minimap2
  samtools

Optional environment variables:
  THREADS        default: 8
  MINIMAP2_BIN   default: minimap2
  SAMTOOLS_BIN   default: samtools

Options:
  --root DIR      permutation root. Default: temp/rrmg_permutation
  --raw           use reads/ rather than qc/
  --no-paf        do not generate read-level PAF.gz
  --force         overwrite existing outputs
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
GENE_REF=""
ROOT="temp/rrmg_permutation"
USE_QC=1
MAKE_PAF=1
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w) WD="$2"; shift 2;;
    -gene|--gene-ref) GENE_REF="$2"; shift 2;;
    --root) ROOT="$2"; shift 2;;
    --raw) USE_QC=0; shift 1;;
    --no-paf) MAKE_PAF=0; shift 1;;
    --force) FORCE=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

: "${WD:?Missing -w}"
: "${GENE_REF:?Missing -gene ref/CYP51A_wt.fasta}"

THREADS="${THREADS:-8}"
MINIMAP2_BIN="${MINIMAP2_BIN:-minimap2}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

WD="$(mkdir -p "$WD" && cd "$WD" && pwd)"
cd "$WD"

[[ -s "$PLAN" ]] || { echo "[ERROR] Plan not found: $PLAN" >&2; exit 1; }
[[ -s "$GENE_REF" ]] || { echo "[ERROR] Gene reference not found: $GENE_REF" >&2; exit 1; }
[[ -s "$FASTQ_MIXER" ]] || { echo "[ERROR] FASTQ audit utility not found: $FASTQ_MIXER" >&2; exit 1; }

marker_matches() {
  local marker="$1" expected="$2"
  [[ -s "$marker" && "$(tr -d '\r\n' < "$marker")" == "$expected" ]]
}

write_marker() {
  local marker="$1" value="$2" temporary="${1}.tmp.$$"
  printf '%s\n' "$value" > "$temporary"
  mv -f -- "$temporary" "$marker"
}

command -v "$MINIMAP2_BIN" >/dev/null 2>&1 || { echo "[ERROR] minimap2 not found: $MINIMAP2_BIN" >&2; exit 1; }
command -v "$SAMTOOLS_BIN" >/dev/null 2>&1 || { echo "[ERROR] samtools not found: $SAMTOOLS_BIN" >&2; exit 1; }

mkdir -p "${WD}/temp"

ROW_ENV="$(mktemp "${WD}/temp/rrmg_readmap_row_env.XXXXXX")"

python3 - "$PLAN" "$ROW_INDEX" "$WD" "$ROOT" "$USE_QC" > "$ROW_ENV" <<'PY'
import csv
import glob
import os
import re
import shlex
import sys
from pathlib import Path

plan = Path(sys.argv[1])
row_index = int(float(sys.argv[2]))
wd = Path(sys.argv[3])
root = sys.argv[4]
use_qc = int(sys.argv[5])

def q(v):
    return shlex.quote(str(v))

def safe(s):
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s).strip())
    return s.strip("_") or "NA"

def resolve(p):
    if not p:
        return ""
    p = Path(p)
    return str(p if p.is_absolute() else wd / p)

root_abs = Path(resolve(root))

chosen = None
with plan.open(newline="", errors="ignore") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)

for row in rows:
    idx = row.get("run_index", "")
    if idx and int(float(idx)) == row_index:
        chosen = row
        break

# Fallback if plan has no run_index: use nth row, 1-based.
if chosen is None and 1 <= row_index <= len(rows):
    chosen = rows[row_index - 1]

if chosen is None:
    raise SystemExit("[ERROR] Row index not found: {}".format(row_index))

seed = chosen.get("seed", "")
sample = chosen.get("sample", "")
run_id = chosen.get("run_id", "")

if not run_id:
    ta = chosen.get("target_abundance", "")
    total_pairs = chosen.get("total_pairs", "")
    if seed and ta and total_pairs:
        run_id = "AirMG_Afum__ta{}__pairs{}__Seed{}".format(ta, total_pairs, seed)
    elif sample and seed:
        run_id = "{}__Seed{}".format(sample, seed)
    else:
        run_id = "row{}".format(row_index)

run_dir = root_abs / run_id

if not run_dir.exists():
    # Robust fallback: glob by seed and abundance/depth.
    seed_pat = "*Seed{}*".format(seed) if seed else "*"
    candidates = list(root_abs.glob(seed_pat))
    if len(candidates) == 1:
        run_dir = candidates[0]
    elif len(candidates) > 1:
        ta = chosen.get("target_abundance", "")
        tp = chosen.get("total_pairs", "")
        candidates2 = [x for x in candidates if ("ta{}".format(ta) in x.name and "pairs{}".format(tp) in x.name)]
        if candidates2:
            run_dir = candidates2[0]

if not run_dir.exists():
    raise SystemExit("[ERROR] Run directory not found: {}".format(run_dir))

read_dir = run_dir / ("qc" if use_qc else "reads")
if not read_dir.exists():
    raise SystemExit("[ERROR] Read directory not found: {}".format(read_dir))

r1s = sorted([
    x for x in read_dir.glob("*_1.fastq.gz")
    if not x.name.endswith(".log") and "target_1" not in x.name
])
if not r1s:
    # fallback: all R1 including target, then choose the longest name? better to fail.
    raise SystemExit("[ERROR] No R1 FASTQ found in {}".format(read_dir))

r1 = r1s[0]
r2 = Path(str(r1).replace("_1.fastq.gz", "_2.fastq.gz"))

if not r2.exists():
    raise SystemExit("[ERROR] Paired R2 FASTQ not found for {}".format(r1))

map_dir = run_dir / "readmap"
map_dir.mkdir(parents=True, exist_ok=True)

total_pairs = chosen.get("total_pairs", "") or chosen.get("computed_total_pairs", "")
target_abundance = chosen.get("target_abundance", "") or "0"
target_pairs = int(round(float(total_pairs) * float(target_abundance))) if total_pairs else ""
background_pairs = int(round(float(total_pairs))) - target_pairs if total_pairs else ""
vals = {
    "ROW_INDEX": row_index,
    "RUN_ID": run_dir.name,
    "RUN_DIR": str(run_dir),
    "READ_DIR": str(read_dir),
    "R1": str(r1),
    "R2": str(r2),
    "MAP_DIR": str(map_dir),
    "SAMPLE": sample,
    "SEED": seed,
    "TARGET_ABUNDANCE": target_abundance,
    "TOTAL_PAIRS": total_pairs,
    "TARGET_PAIRS_EXPECTED": target_pairs,
    "BACKGROUND_PAIRS_EXPECTED": background_pairs,
}
for k, v in vals.items():
    print("{}={}".format(k, q(v)))
PY

# shellcheck disable=SC1090
source "$ROW_ENV"
rm -f "$ROW_ENV"

PREFIX="${MAP_DIR}/reads_to_GENE"
BAM="${PREFIX}.bam"
DEPTH="${PREFIX}.depth.tsv"
SUMMARY="${PREFIX}.coverage_summary.tsv"
PAF_GZ="${PREFIX}.paf.gz"
OUTPUT_PROVENANCE="${PREFIX}.mix_generation.txt"

mkdir -p "$MAP_DIR"

echo "[INFO] ROW_INDEX         = $ROW_INDEX"
echo "[INFO] RUN_ID            = $RUN_ID"
echo "[INFO] RUN_DIR           = $RUN_DIR"
echo "[INFO] R1                = $R1"
echo "[INFO] R2                = $R2"
echo "[INFO] GENE_REF          = $GENE_REF"
echo "[INFO] THREADS           = $THREADS"
echo "[INFO] MINIMAP2_BIN      = $MINIMAP2_BIN"
echo "[INFO] SAMTOOLS_BIN      = $SAMTOOLS_BIN"
echo "[INFO] MAKE_PAF          = $MAKE_PAF"

[[ -n "$TOTAL_PAIRS" ]] || { echo "[ERROR] Plan row has no total_pairs/computed_total_pairs" >&2; exit 1; }
MIX_AUDIT="${RUN_DIR}/design/read_mix_audit.tsv"
if [[ "$USE_QC" -eq 1 ]]; then
  python3 "$FASTQ_MIXER" validate-mix --audit "$MIX_AUDIT" --expected-pairs "$TOTAL_PAIRS" \
    --expected-background-pairs "$BACKGROUND_PAIRS_EXPECTED" --expected-target-pairs "$TARGET_PAIRS_EXPECTED"
  MIX_GENERATION="$(awk -F '\t' '$1=="mix_generation_id"{print $2}' "$MIX_AUDIT")"
  FASTP_JSON="${R1%_1.fastq.gz}.fastp.json"
  QC_PROVENANCE="${R1%_1.fastq.gz}.mix_generation.txt"
  python3 "$FASTQ_MIXER" validate-fastp --json "$FASTP_JSON" --expected-pairs "$TOTAL_PAIRS"
  marker_matches "$QC_PROVENANCE" "$MIX_GENERATION" || {
      echo "[ERROR] QC reads/fastp report belong to an older/unaudited Step 2 mix." >&2
      echo "        Rerun corrected Step 3 before supplementary read mapping." >&2
      exit 1
    }
else
  python3 "$FASTQ_MIXER" validate-mix \
    --audit "$MIX_AUDIT" --expected-pairs "$TOTAL_PAIRS" \
    --expected-background-pairs "$BACKGROUND_PAIRS_EXPECTED" --expected-target-pairs "$TARGET_PAIRS_EXPECTED" \
    --r1 "$R1" --r2 "$R2"
  MIX_GENERATION="$(awk -F '\t' '$1=="mix_generation_id"{print $2}' "$MIX_AUDIT")"
fi

if [[ "$FORCE" -eq 0 && -s "$SUMMARY" && -s "$DEPTH" ]]; then
  OUTPUTS_CURRENT=1
  marker_matches "$OUTPUT_PROVENANCE" "$MIX_GENERATION" || OUTPUTS_CURRENT=0
  if [[ "$MAKE_PAF" -eq 1 ]]; then
    [[ -s "$PAF_GZ" ]] || OUTPUTS_CURRENT=0
  fi
  if [[ "$OUTPUTS_CURRENT" -eq 1 ]]; then
    echo "[SKIP] existing audited read-depth outputs: $SUMMARY"
    exit 0
  fi
  echo "[STALE] Existing Step 4 outputs do not match the current read-mix generation; regenerating them."
fi

rm -f "$BAM" "${BAM}.bai" "$DEPTH" "$SUMMARY" "$PAF_GZ" "$OUTPUT_PROVENANCE"

echo "[MAP] reads -> target gene BAM"

"$MINIMAP2_BIN" \
  -t "$THREADS" \
  -ax sr \
  --secondary=no \
  "$GENE_REF" "$R1" "$R2" \
  2> "${PREFIX}.bam.minimap2.log" \
| "$SAMTOOLS_BIN" view -@ "$THREADS" -b -F 4 - \
| "$SAMTOOLS_BIN" sort -@ "$THREADS" -o "$BAM" -

"$SAMTOOLS_BIN" index "$BAM"

echo "[DEPTH] samtools depth -aa"
"$SAMTOOLS_BIN" depth -aa "$BAM" > "$DEPTH"

if [[ "$MAKE_PAF" -eq 1 ]]; then
  echo "[PAF] reads -> target gene PAF.gz"
  "$MINIMAP2_BIN" \
    -t "$THREADS" \
    -x sr \
    --secondary=no \
    "$GENE_REF" "$R1" "$R2" \
    2> "${PREFIX}.paf.minimap2.log" \
  | gzip -c > "$PAF_GZ"
fi

MAPPED_READS="$("$SAMTOOLS_BIN" view -c "$BAM" || echo 0)"

python3 - "$DEPTH" "$SUMMARY" "$GENE_REF" "$RUN_ID" "$ROW_INDEX" "$SEED" "$TARGET_ABUNDANCE" "$TOTAL_PAIRS" "$R1" "$R2" "$BAM" "$PAF_GZ" "$MAPPED_READS" <<'PY'
import csv
import gzip
import statistics
import sys
from pathlib import Path

depth_file, out_file, gene_ref, run_id, row_index, seed, ta, total_pairs, r1, r2, bam, paf_gz, mapped_reads = sys.argv[1:]

def open_text(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "r", errors="ignore")

def fasta_len(path):
    n = 0
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n

depths = []
target_names = set()
with open(depth_file, errors="ignore") as fh:
    for line in fh:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        target_names.add(parts[0])
        depths.append(int(float(parts[2])))

gene_len = fasta_len(gene_ref)
if not depths and gene_len > 0:
    depths = [0] * gene_len

covered_bp = sum(1 for x in depths if x > 0)
ge5_bp = sum(1 for x in depths if x >= 5)
ge10_bp = sum(1 for x in depths if x >= 10)
mean_depth = sum(depths) / float(len(depths)) if depths else 0
median_depth = statistics.median(depths) if depths else 0
max_depth = max(depths) if depths else 0
breadth_pct = covered_bp / float(len(depths)) * 100 if depths else 0
ge5_pct = ge5_bp / float(len(depths)) * 100 if depths else 0
ge10_pct = ge10_bp / float(len(depths)) * 100 if depths else 0

fields = [
    "RunID", "RowIndex", "Seed", "TargetAbundanceDesign", "TotalPairs",
    "TargetGeneFA", "TargetGeneBp", "MappedReadAlignments",
    "CoveredBp", "BreadthPct", "CoveredBpDepthGE5", "BreadthPctDepthGE5",
    "CoveredBpDepthGE10", "BreadthPctDepthGE10",
    "MeanDepth", "MedianDepth", "MaxDepth",
    "DepthTSV", "BAM", "PAFGZ", "R1", "R2"
]

row = {
    "RunID": run_id,
    "RowIndex": row_index,
    "Seed": seed,
    "TargetAbundanceDesign": ta,
    "TotalPairs": total_pairs,
    "TargetGeneFA": gene_ref,
    "TargetGeneBp": gene_len,
    "MappedReadAlignments": mapped_reads,
    "CoveredBp": covered_bp,
    "BreadthPct": "{:.6f}".format(breadth_pct),
    "CoveredBpDepthGE5": ge5_bp,
    "BreadthPctDepthGE5": "{:.6f}".format(ge5_pct),
    "CoveredBpDepthGE10": ge10_bp,
    "BreadthPctDepthGE10": "{:.6f}".format(ge10_pct),
    "MeanDepth": "{:.6f}".format(mean_depth),
    "MedianDepth": "{:.6f}".format(median_depth),
    "MaxDepth": max_depth,
    "DepthTSV": depth_file,
    "BAM": bam,
    "PAFGZ": paf_gz if Path(paf_gz).exists() else "",
    "R1": r1,
    "R2": r2,
}

with open(out_file, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerow(row)
PY

[[ -s "$SUMMARY" ]] || { echo "[ERROR] Summary not created: $SUMMARY" >&2; exit 1; }
write_marker "$OUTPUT_PROVENANCE" "$MIX_GENERATION"
echo "[DONE] $SUMMARY"
