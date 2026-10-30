#!/usr/bin/env bash
# Add evidence from existing contigs/BAMs; no FASTQ, QC, assembly or read remapping.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$SCRIPT_DIR/rrmg_gene_evidence.py" backfill "$@"
