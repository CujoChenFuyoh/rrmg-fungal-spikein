#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${RRMG_CONFIG:-$SCRIPT_DIR/rrmg_config.sh}"
[[ -r "$CONFIG" ]] || { echo "[ERROR] Configuration not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"
: "${PIPELINE_DIR:?Set PIPELINE_DIR}"
: "${WD:?Set WD}"
source "$PIPELINE_DIR/rrmg_environment.sh"
errors=0
check() {
  if command -v "$1" >/dev/null 2>&1; then
    printf '[OK] %s: %s\n' "$1" "$(command -v "$1")"
  else
    echo "[ERROR] Missing tool in active PATH: $1" >&2
    errors=$((errors + 1))
  fi
}
for tool in python3 datasets unzip java fastp minimap2 samtools metaspades.py; do check "$tool"; done
[[ "${RUN_KRAKEN:-0}" == 1 ]] && check kraken2
if ! command -v "${VAMB_COMMAND:-vamb}" >/dev/null 2>&1; then
  echo '[WARN] Optional VAMB missing; required only for step3b.'
fi
if [[ -s "${RANDOMREADSMG_SH:-}" ]]; then
  echo "[OK] RandomReadsMG: $RANDOMREADSMG_SH"
else
  echo '[ERROR] RandomReadsMG missing: put randomreadsmg.sh on PATH or set RANDOMREADSMG_SH.' >&2
  errors=$((errors + 1))
fi
if [[ -s "${RANDOMREADS_SH:-}" ]]; then
  echo "[OK] RandomReads: $RANDOMREADS_SH"
else
  echo '[ERROR] RandomReads missing: put randomreads.sh on PATH or set RANDOMREADS_SH.' >&2
  errors=$((errors + 1))
fi
[[ -z "${JAVA_BIN:-}" ]] || echo "[INFO] Java real path: $(readlink -f "$JAVA_BIN")"
for file in rrmg_submit.sh rrmg_run_local.sh rrmg_run_one.sh rrmg_environment.sh \
  rrmg_step1_make_background_reads.sh rrmg_step2_make_target_and_merge.sh \
  rrmg_mix_paired_fastq.py \
  rrmg_step3_downstream_analysis.sh rrmg_step4_5_supplementary_downstream.sh \
  rrmg_step4_reads_to_gene_depth.sh rrmg_step5_afum_target_reads_genome_cov.sh; do
  if [[ ! -s "$PIPELINE_DIR/$file" ]]; then
    echo "[ERROR] Missing core file: $file" >&2
    errors=$((errors + 1))
  fi
done
for variable in PLAN BACKGROUND_PLAN BACKGROUND_INDEX GENE_REF TARGET_GENOME; do
  value="${!variable:-}"
  [[ "$value" == /* ]] || value="$WD/$value"
  [[ -s "$value" ]] || echo "[WARN] $variable not present yet: $value"
done
echo "[INFO] No environment was activated, no scheduler called, no analysis run."
echo "[SUMMARY] errors=$errors (missing future inputs are warnings)"
[[ "$errors" == 0 ]]
