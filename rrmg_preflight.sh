#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${RRMG_CONFIG:-${SCRIPT_DIR}/rrmg_config.sh}"

if [[ ! -r "$CONFIG" ]]; then
  echo "[ERROR] Configuration not found: $CONFIG" >&2
  echo "        cp ${SCRIPT_DIR}/rrmg_config.example.sh ${SCRIPT_DIR}/rrmg_config.sh" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG"

: "${PIPELINE_DIR:?Set PIPELINE_DIR in $CONFIG}"
: "${WD:?Set WD in $CONFIG}"
: "${CONDA_SH:?Set CONDA_SH in $CONFIG}"
: "${RRMG_CONDA_ENV:?Set RRMG_CONDA_ENV in $CONFIG}"
: "${RANDOMREADSMG_SH:?Set RANDOMREADSMG_SH in $CONFIG}"
: "${JAVA_BIN:?Set JAVA_BIN in $CONFIG}"

errors=0
warnings=0

ok() { printf '[OK] %s\n' "$*"; }
warn() { printf '[WARN] %s\n' "$*" >&2; warnings=$((warnings + 1)); }
fail() { printf '[ERROR] %s\n' "$*" >&2; errors=$((errors + 1)); }

[[ -d "$PIPELINE_DIR" ]] && ok "pipeline directory: $PIPELINE_DIR" || fail "missing PIPELINE_DIR: $PIPELINE_DIR"
[[ -d "$WD" ]] && ok "working directory: $WD" || fail "missing WD: $WD"
[[ -r "$CONDA_SH" ]] && ok "conda initialization: $CONDA_SH" || fail "missing CONDA_SH: $CONDA_SH"
[[ -s "$RANDOMREADSMG_SH" ]] && ok "RandomReadsMG: $RANDOMREADSMG_SH" || fail "missing RANDOMREADSMG_SH: $RANDOMREADSMG_SH"
[[ -x "$JAVA_BIN" ]] && ok "Java: $JAVA_BIN" || fail "JAVA_BIN is not executable: $JAVA_BIN"

required_scripts=(
  rrmg_run_local.sh
  rrmg_submit.sh
  rrmg_make_plan.py
  rrmg_make_background_plan.py
  rrmg_step1_make_background_reads.sh
  rrmg_step2_make_target_and_merge.sh
  rrmg_step3_downstream_analysis.sh
  rrmg_step3b_vamb_binning.sh
  rrmg_evaluate_vamb_bins.py
  rrmg_step1_background_array.pbs
  rrmg_step2_target_merge_array.pbs
  rrmg_step3_downstream_array.pbs
  rrmg_step3b_vamb_array.pbs
  rrmg_step4_5_supplementary_downstream_array.pbs
)

for script in "${required_scripts[@]}"; do
  [[ -s "$PIPELINE_DIR/$script" ]] || fail "missing pipeline script: $PIPELINE_DIR/$script"
done

if [[ -r "$CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH"
  if conda run -n "$RRMG_CONDA_ENV" python3 --version >/dev/null 2>&1; then
    ok "Conda environment: $RRMG_CONDA_ENV"
  else
    fail "cannot run Python in Conda environment: $RRMG_CONDA_ENV"
  fi

  tools=(python3 fastp minimap2 samtools metaspades.py vamb)
  [[ "${RUN_KRAKEN:-0}" -eq 1 ]] && tools+=(kraken2)
  for tool in "${tools[@]}"; do
    if conda run -n "$RRMG_CONDA_ENV" bash -lc "command -v '$tool'" >/dev/null 2>&1; then
      ok "tool in ${RRMG_CONDA_ENV}: $tool"
    else
      fail "tool missing from ${RRMG_CONDA_ENV}: $tool"
    fi
  done
fi

resolve_from_wd() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s/%s\n' "$WD" "$1" ;;
  esac
}

for variable in PLAN BACKGROUND_PLAN BACKGROUND_INDEX GENE_REF TARGET_GENOME; do
  value="${!variable:-}"
  [[ -n "$value" ]] || { warn "$variable is empty"; continue; }
  path="$(resolve_from_wd "$value")"
  [[ -s "$path" ]] && ok "$variable: $path" || warn "$variable not present yet: $path"
done

if [[ "${RUN_KRAKEN:-0}" -eq 1 ]]; then
  [[ -d "${KRAKEN_DB:-}" ]] && ok "Kraken2 database: $KRAKEN_DB" || fail "RUN_KRAKEN=1 but KRAKEN_DB is missing: ${KRAKEN_DB:-unset}"
fi

printf '\nPreflight summary: errors=%d warnings=%d\n' "$errors" "$warnings"
[[ "$errors" -eq 0 ]]
