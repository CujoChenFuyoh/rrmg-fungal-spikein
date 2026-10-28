#!/usr/bin/env bash
# Internal single-row dispatcher. No scheduler calls or environment activation.
set -euo pipefail
STEP="${1:?Specify step}"
ROW="${2:?Specify original plan index}"
[[ "$ROW" =~ ^[1-9][0-9]*$ ]] || { echo '[ERROR] Invalid row index' >&2; exit 1; }
source "${RRMG_CONFIG:?Set RRMG_CONFIG}"
: "${WD:?Set WD in config}"
: "${PIPELINE_DIR:?Set PIPELINE_DIR in config}"
source "$PIPELINE_DIR/rrmg_environment.sh"
cd "$WD"
mkdir -p temp result
PLAN="${PLAN:-input/rrmg/rrmg_permutation_plan.tsv}"
GENE_REF="${GENE_REF:-ref/CYP51A_wt.fasta}"
TARGET_GENOME="${TARGET_GENOME:-ref/Afum_Af293.fasta}"
ARGS=()
force_arg() { [[ "${RRMG_FORCE_OVERRIDE:-$1}" == 1 ]] && ARGS+=(--force); return 0; }
require_file() { [[ -s "$1" ]] || { echo "[ERROR] Missing input: $1" >&2; exit 1; }; }
if [[ -n "${RRMG_THREADS_OVERRIDE:-}" && ! "${RRMG_THREADS_OVERRIDE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] RRMG_THREADS_OVERRIDE must be a positive integer" >&2
  exit 1
fi
if [[ -n "${RRMG_MEMORY_GB_OVERRIDE:-}" && ! "${RRMG_MEMORY_GB_OVERRIDE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] RRMG_MEMORY_GB_OVERRIDE must be a positive integer in GiB" >&2
  exit 1
fi
if [[ -n "${RRMG_MEMORY_GB_OVERRIDE:-}" ]]; then
  case "$STEP" in
    step1|step2|step3) ;;
    *) echo "[ERROR] Memory override is unsupported for $STEP; set scheduler memory instead" >&2; exit 1;;
  esac
fi
case "$STEP" in
  step1)
    require_file "${RANDOMREADSMG_SH:-randomreadsmg.sh not in PATH; configure RANDOMREADSMG_SH}"
    command -v "${JAVA_BIN:-java}" >/dev/null || { echo '[ERROR] Java missing from active environment' >&2; exit 1; }
    ;;
  step2)
    require_file "${RANDOMREADS_SH:-randomreads.sh not in PATH; configure RANDOMREADS_SH}"
    command -v "${JAVA_BIN:-java}" >/dev/null || { echo '[ERROR] Java missing from active environment' >&2; exit 1; }
    ;;
esac
case "$STEP" in
  step1)
    SCRIPT=rrmg_step1_make_background_reads.sh
    export THREADS="${RRMG_THREADS_OVERRIDE:-${STEP1_THREADS:-4}}"
    export BBMAP_XMX="${RRMG_MEMORY_GB_OVERRIDE:+${RRMG_MEMORY_GB_OVERRIDE}g}"
    BBMAP_XMX="${BBMAP_XMX:-${STEP1_BBMAP_XMX:-4g}}"
    export BBMAP_XMX BACKGROUND_SEED="${BACKGROUND_SEED:-1000}"
    ARGS=(--background-plan "${BACKGROUND_PLAN:-temp/rrmg_background/background_plan.tsv}" --row-index "$ROW" -w "$WD")
    force_arg "${STEP1_FORCE:-0}"
    ;;
  step2)
    SCRIPT=rrmg_step2_make_target_and_merge.sh
    export THREADS="${RRMG_THREADS_OVERRIDE:-${STEP2_THREADS:-1}}"
    export BBMAP_XMX="${RRMG_MEMORY_GB_OVERRIDE:+${RRMG_MEMORY_GB_OVERRIDE}g}"
    BBMAP_XMX="${BBMAP_XMX:-${STEP2_BBMAP_XMX:-4g}}"
    export BBMAP_XMX
    ARGS=("$PLAN" "$ROW" -w "$WD" --background-index "${BACKGROUND_INDEX:-temp/rrmg_background/background_index.tsv}")
    [[ "${STEP2_KEEP_TARGET_READS:-1}" == 1 ]] && ARGS+=(--keep-target-reads)
    force_arg "${STEP2_FORCE:-0}"
    ;;
  step3)
    SCRIPT=rrmg_step3_downstream_analysis.sh
    require_file "$GENE_REF"
    ARGS=("$PLAN" "$ROW" -w "$WD" -gene "$GENE_REF" -t "${RRMG_THREADS_OVERRIDE:-${STEP3_THREADS:-4}}" -mem "${RRMG_MEMORY_GB_OVERRIDE:-${STEP3_SPADES_MEM_GB:-16}}" --phred-offset "${STEP3_PHRED_OFFSET:-33}")
    case "${STEP3_ONLY_ASSEMBLER:-1}" in
      1) ARGS+=(--only-assembler);;
      0) ARGS+=(--with-error-correction);;
      *) echo '[ERROR] STEP3_ONLY_ASSEMBLER must be 0 or 1' >&2; exit 1;;
    esac
    if [[ "${RUN_KRAKEN:-0}" == 1 ]]; then
      ARGS+=(--db "${KRAKEN_DB:?Set KRAKEN_DB}")
    else
      ARGS+=(--skip-kraken)
    fi
    [[ "${STEP3_SKIP_QC:-0}" == 1 ]] && ARGS+=(--skip-qc)
    [[ "${STEP3_SKIP_ASSEMBLY:-0}" == 1 ]] && ARGS+=(--skip-assemble)
    force_arg "${STEP3_FORCE:-0}"
    ;;
  step3b)
    SCRIPT=rrmg_step3b_vamb_binning.sh
    ARGS=("$PLAN" "$ROW" -w "$WD" -t "${RRMG_THREADS_OVERRIDE:-${STEP3B_THREADS:-4}}" --min-contig "${VAMB_MIN_CONTIG:-2000}" --min-bin-size "${VAMB_MIN_BIN_SIZE:-200000}" --vamb-seed "${VAMB_SEED:-42}" --vamb-command "${VAMB_COMMAND:-vamb}")
    [[ "${VAMB_USE_CUDA:-0}" == 1 ]] && ARGS+=(--cuda)
    force_arg "${VAMB_FORCE:-0}"
    ;;
  step4|step5|step45)
    SCRIPT=rrmg_step4_5_supplementary_downstream.sh
    export THREADS="${RRMG_THREADS_OVERRIDE:-${STEP45_THREADS:-4}}" MINIMAP2_BIN="${MINIMAP2_BIN:-minimap2}" SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
    export READS_TO_GENE_SCRIPT="$PIPELINE_DIR/rrmg_step4_reads_to_gene_depth.sh"
    export TARGET_GENOME_COV_SCRIPT="$PIPELINE_DIR/rrmg_step5_afum_target_reads_genome_cov.sh"
    ARGS=("$PLAN" "$ROW" -w "$WD" -gene "$GENE_REF" --target-genome "$TARGET_GENOME")
    [[ "$STEP" == step4 ]] && ARGS+=(--skip-target-genome)
    [[ "$STEP" == step5 ]] && ARGS+=(--skip-read-to-gene)
    [[ "${STEP4_MAP_TARGET_ONLY:-1}" == 0 && "$STEP" != step5 ]] && ARGS+=(--skip-target-only-gene)
    [[ "${STEP45_MAKE_BINS:-1}" == 1 ]] && ARGS+=(--make-bins --bin-size "${STEP45_BIN_SIZE:-10000}")
    force_arg "${STEP45_FORCE:-0}"
    ;;
  *) echo "[ERROR] Unknown step: $STEP" >&2; exit 1;;
esac
printf '[COMMAND]'; printf ' %q' bash "$PIPELINE_DIR/$SCRIPT" "${ARGS[@]}"; printf '\n'
exec bash "$PIPELINE_DIR/$SCRIPT" "${ARGS[@]}"
