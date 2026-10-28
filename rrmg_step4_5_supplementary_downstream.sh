#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  rrmg_step4_5_supplementary_downstream.sh plan.tsv ROW_INDEX \
    -w /path/to/project \
    -gene ref/CYP51A_wt.fasta \
    --target-genome ref/Afum_Af293.fasta

Purpose:
  Run both supplementary downstream checks for one RRMG permutation row:
    A. spiked metagenome and target-only reads -> target gene, with per-base depth
    B. target-only simulated Afum reads -> Afum genome, with genome-wide coverage

This wrapper calls the existing working scripts:
  rrmg_step4_reads_to_gene_depth.sh
  rrmg_step5_afum_target_reads_genome_cov.sh

Options:
  -w DIR                         working directory
  -gene FILE, --gene-ref FILE    target gene FASTA, e.g. CYP51A_wt.fasta
  --target-genome FILE           target genome FASTA, e.g. Afum_Af293.fasta
  --root DIR                     permutation root, default temp/rrmg_permutation

  --skip-read-to-gene            skip spiked metagenome reads -> target gene
  --skip-target-genome           skip target-only reads -> Afum genome

  --raw                          Step4: use raw merged reads instead of QC reads
  --skip-target-only-gene        Step4: skip target-only reads -> target gene
  --no-paf                       Step4: skip read-level PAF.gz
  --make-bins                    Step5: generate binned genome coverage table
  --bin-size N                   Step5: bin size, default 10000
  --remove-target-genome-bam     Step5: remove genome BAM/BAI after summary
  --force                        overwrite existing Step4 and Step5 outputs
  --readlen N                    Step5 fallback read length, default 150

Environment:
  THREADS                        minimap2/samtools threads, default 8
  MINIMAP2_BIN                   default minimap2
  SAMTOOLS_BIN                   default samtools
  READS_TO_GENE_SCRIPT           default rrmg_step4_reads_to_gene_depth.sh
  TARGET_GENOME_COV_SCRIPT       default rrmg_step5_afum_target_reads_genome_cov.sh
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
TARGET_GENOME=""
ROOT="temp/rrmg_permutation"

SKIP_READ_TO_GENE=0
SKIP_TARGET_GENOME=0

USE_RAW=0
SKIP_TARGET_ONLY_GENE=0
NO_PAF=0
MAKE_BINS=0
BIN_SIZE=10000
REMOVE_TARGET_GENOME_BAM=0
FORCE=0
READLEN_DEFAULT=150

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w) WD="$2"; shift 2;;
    -gene|--gene-ref) GENE_REF="$2"; shift 2;;
    --target-genome) TARGET_GENOME="$2"; shift 2;;
    --root) ROOT="$2"; shift 2;;
    --skip-read-to-gene) SKIP_READ_TO_GENE=1; shift 1;;
    --skip-target-genome) SKIP_TARGET_GENOME=1; shift 1;;
    --raw) USE_RAW=1; shift 1;;
    --skip-target-only-gene) SKIP_TARGET_ONLY_GENE=1; shift 1;;
    --no-paf) NO_PAF=1; shift 1;;
    --make-bins) MAKE_BINS=1; shift 1;;
    --bin-size) BIN_SIZE="$2"; shift 2;;
    --remove-target-genome-bam) REMOVE_TARGET_GENOME_BAM=1; shift 1;;
    --force) FORCE=1; shift 1;;
    --readlen) READLEN_DEFAULT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

: "${WD:?Missing -w}"

if [[ "$SKIP_READ_TO_GENE" -eq 0 ]]; then
  : "${GENE_REF:?Missing -gene ref/CYP51A_wt.fasta}"
fi

THREADS="${THREADS:-8}"
MINIMAP2_BIN="${MINIMAP2_BIN:-minimap2}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

READS_TO_GENE_SCRIPT="${READS_TO_GENE_SCRIPT:-rrmg_step4_reads_to_gene_depth.sh}"
TARGET_GENOME_COV_SCRIPT="${TARGET_GENOME_COV_SCRIPT:-rrmg_step5_afum_target_reads_genome_cov.sh}"

WD="$(mkdir -p "$WD" && cd "$WD" && pwd)"
cd "$WD"

[[ -s "$PLAN" ]] || { echo "[ERROR] Missing PLAN: $PLAN" >&2; exit 1; }

command -v "$MINIMAP2_BIN" >/dev/null 2>&1 || { echo "[ERROR] minimap2 not found: $MINIMAP2_BIN" >&2; exit 1; }
command -v "$SAMTOOLS_BIN" >/dev/null 2>&1 || { echo "[ERROR] samtools not found: $SAMTOOLS_BIN" >&2; exit 1; }

echo "++++++++++++++++++++++++++++++++++++++++"
echo "Supplementary downstream analysis"
echo "ROW_INDEX              = $ROW_INDEX"
echo "WD                     = $WD"
echo "PLAN                   = $PLAN"
echo "GENE_REF               = ${GENE_REF:-NA}"
echo "TARGET_GENOME          = ${TARGET_GENOME:-from plan/design}"
echo "ROOT                   = $ROOT"
echo "THREADS                = $THREADS"
echo "SKIP_READ_TO_GENE      = $SKIP_READ_TO_GENE"
echo "SKIP_TARGET_GENOME     = $SKIP_TARGET_GENOME"
echo "SKIP_TARGET_ONLY_GENE  = $SKIP_TARGET_ONLY_GENE"
echo "++++++++++++++++++++++++++++++++++++++++"

if [[ "$SKIP_READ_TO_GENE" -eq 0 ]]; then
  [[ -s "$READS_TO_GENE_SCRIPT" ]] || { echo "[ERROR] Missing script: $READS_TO_GENE_SCRIPT" >&2; exit 1; }
  [[ -s "$GENE_REF" ]] || { echo "[ERROR] Missing gene reference: $GENE_REF" >&2; exit 1; }

  CMD=(bash "$READS_TO_GENE_SCRIPT" "$PLAN" "$ROW_INDEX" -w "$WD" -gene "$GENE_REF" --root "$ROOT")
  [[ "$USE_RAW" -eq 1 ]] && CMD+=(--raw)
  [[ "$SKIP_TARGET_ONLY_GENE" -eq 1 ]] && CMD+=(--skip-target-only)
  [[ "$NO_PAF" -eq 1 ]] && CMD+=(--no-paf)
  [[ "$FORCE" -eq 1 ]] && CMD+=(--force)

  echo "[RUN-A] spiked metagenome + target-only reads -> target gene"
  printf ' %q' "${CMD[@]}"; echo
  "${CMD[@]}"
else
  echo "[SKIP-A] read-to-gene mapping skipped"
fi

if [[ "$SKIP_TARGET_GENOME" -eq 0 ]]; then
  [[ -s "$TARGET_GENOME_COV_SCRIPT" ]] || { echo "[ERROR] Missing script: $TARGET_GENOME_COV_SCRIPT" >&2; exit 1; }

  CMD=(bash "$TARGET_GENOME_COV_SCRIPT" "$PLAN" "$ROW_INDEX" -w "$WD" --root "$ROOT" --readlen "$READLEN_DEFAULT")
  [[ -n "$TARGET_GENOME" ]] && CMD+=(--target-genome "$TARGET_GENOME")
  [[ "$MAKE_BINS" -eq 1 ]] && CMD+=(--make-bins --bin-size "$BIN_SIZE")
  [[ "$REMOVE_TARGET_GENOME_BAM" -eq 1 ]] && CMD+=(--remove-bam)
  [[ "$FORCE" -eq 1 ]] && CMD+=(--force)

  echo "[RUN-B] target-only Afum reads -> Afum genome"
  printf ' %q' "${CMD[@]}"; echo
  "${CMD[@]}"
else
  echo "[SKIP-B] target-only genome coverage skipped"
fi

echo "++++++++++++++++++++++++++++++++++++++++"
echo "[DONE] Supplementary downstream ROW_INDEX=$ROW_INDEX"
echo "++++++++++++++++++++++++++++++++++++++++"
