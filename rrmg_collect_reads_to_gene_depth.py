#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Collect all per-run CYP51A read-depth summaries into one TSV.

Collect the available ``*/readmap/reads_to_GENE.coverage_summary.tsv`` files.
This merges existing summaries rather than calculating coverage. Partial cohorts
are allowed; compare collected RunIDs with the plan before interpreting results.
"""

import argparse
import csv
import os
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Collect RRMG reads-to-gene depth summaries.")
    p.add_argument("--root", default="temp/rrmg_permutation", help="Permutation root directory.")
    p.add_argument("--out", default="result/rrmg_reads_to_gene_depth_summary.tsv", help="Output TSV.")
    return p.parse_args()


def main():
    args = parse_args()
    root = Path(args.root)
    files = sorted(root.glob("*/readmap/reads_to_GENE.coverage_summary.tsv"))
    if not files:
        raise SystemExit("[ERROR] No readmap coverage summaries found under {}".format(root))

    rows = []
    fields = None

    for f in files:
        with f.open(newline="", errors="ignore") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if fields is None:
                    fields = reader.fieldnames
                rows.append(row)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print("[OK] collected {} summaries".format(len(rows)))
    print("[OK] wrote {}".format(out))


if __name__ == "__main__":
    main()
