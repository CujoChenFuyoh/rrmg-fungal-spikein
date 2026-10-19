#!/usr/bin/env bash
set -euo pipefail

# Final collector for the optional combined Step 4+5 analysis.
# Positional arguments:
#   1) project working directory
#   2) permutation run root, relative to the project unless absolute
# Collect existing summaries; audit missing IDs against the plan for partial cohorts.

WD="${1:-$PWD}"
ROOT="${2:-temp/rrmg_permutation}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$WD"
mkdir -p result

[[ -s "$SCRIPT_DIR/rrmg_collect_reads_to_gene_depth.py" ]] || { echo "[ERROR] Missing rrmg_collect_reads_to_gene_depth.py"; exit 1; }
[[ -s "$SCRIPT_DIR/rrmg_collect_afum_genome_cov.py" ]] || { echo "[ERROR] Missing rrmg_collect_afum_genome_cov.py"; exit 1; }

python3 "$SCRIPT_DIR/rrmg_collect_reads_to_gene_depth.py" \
  --root "$ROOT" \
  --out result/rrmg_reads_to_gene_depth_summary.tsv

python3 "$SCRIPT_DIR/rrmg_collect_afum_genome_cov.py" \
  --root "$ROOT" \
  --out result/rrmg_afum_target_reads_genome_coverage_summary.tsv

echo "[DONE] Supplementary downstream summaries:"
echo "  result/rrmg_reads_to_gene_depth_summary.tsv"
echo "  result/rrmg_afum_target_reads_genome_coverage_summary.tsv"
