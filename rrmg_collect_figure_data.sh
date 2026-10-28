#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'USAGE'
Usage:
  rrmg_collect_figure_data.sh --wd WORKDIR --gene GENE_FASTA [options]

Collect all completed Step3-5 plotting inputs into WORKDIR/figure_data.
No mapping, read simulation, assembly or R plotting is performed.

Options:
  --wd DIR             RRMG working directory (required)
  --plan FILE          plan path relative to WD or absolute
                       default: input/rrmg/rrmg_permutation_plan.tsv
  --root DIR           permutation root relative to WD or absolute
                       default: temp/rrmg_permutation
  --gene FILE          target-gene FASTA (required)
  --out DIR            managed output directory; default: WD/figure_data
  --sites FILE         optional validated-site source TSV
  --catalogue FILE     optional catalogue TSV; required with --sites
USAGE
}

WD=""
PLAN="input/rrmg/rrmg_permutation_plan.tsv"
ROOT="temp/rrmg_permutation"
GENE=""
OUT=""
SITES=""
CATALOGUE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --wd) WD="$2"; shift 2;;
    --plan) PLAN="$2"; shift 2;;
    --root) ROOT="$2"; shift 2;;
    --gene) GENE="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --sites) SITES="$2"; shift 2;;
    --catalogue) CATALOGUE="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage; exit 1;;
  esac
done

: "${WD:?Use --wd WORKDIR}"
: "${GENE:?Use --gene GENE_FASTA}"
if [[ -n "$SITES" || -n "$CATALOGUE" ]]; then
  [[ -n "$SITES" && -n "$CATALOGUE" ]] || { echo '[ERROR] Use --sites and --catalogue together' >&2; exit 1; }
fi

WD="$(cd "$WD" && pwd)"
cd "$WD"
[[ -s "$GENE" ]] || { echo "[ERROR] Missing gene FASTA: $GENE" >&2; exit 1; }
[[ -s "$PLAN" ]] || { echo "[ERROR] Missing plan: $PLAN" >&2; exit 1; }
OUT="${OUT:-$WD/figure_data}"
mkdir -p result

python3 "$SCRIPT_DIR/rrmg_collect_map_results.py" \
  --summary "$WD/result/rrmg_split_pipeline_summary.tsv" \
  --plan "$PLAN" \
  --root "$ROOT" \
  --out "$WD/result/rrmg_target_gene_map_merged.tsv" \
  --wd "$WD"

bash "$SCRIPT_DIR/rrmg_collect_supplementary_downstream.sh" "$WD" "$ROOT"

CMD=(python3 "$SCRIPT_DIR/rrmg_collect_figure_data.py"
  --wd "$WD" --plan "$PLAN" --root "$ROOT" --gene "$GENE" --out "$OUT")
if [[ -n "$SITES" ]]; then
  CMD+=(--sites "$SITES" --catalogue "$CATALOGUE")
fi
printf '[COMMAND]'; printf ' %q' "${CMD[@]}"; printf '\n'
"${CMD[@]}"

test -s "$OUT/COMPLETE.json"
test -s "$OUT/FILE_MANIFEST.tsv"
echo "[DONE] Download this complete directory to the local Windows work directory: $OUT"
