#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Collect all target-only Afum genome-coverage summaries into one TSV.

Run this after Step 5 (or the combined Step 4+5 wrapper) has completed. This is
a collector only: the mapping and coverage calculations are performed by
``rrmg_step5_afum_target_reads_genome_cov.sh``.
"""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect simulated Afum target-only reads to genome coverage summaries."
    )
    parser.add_argument("--root", default="temp/rrmg_permutation")
    parser.add_argument(
        "--out",
        default="result/rrmg_afum_target_reads_genome_coverage_summary.tsv",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    files = sorted(
        Path(args.root).glob(
            "*/target_readmap/afum_target_reads_to_genome.summary.tsv"
        )
    )
    if not files:
        raise SystemExit(
            "[ERROR] No target_readmap summaries found under {}".format(args.root)
        )

    rows = []
    fields = None
    for path in files:
        with path.open(newline="", errors="ignore") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if fields is None:
                    fields = reader.fieldnames
                rows.append(row)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    print("[OK] collected {} summaries".format(len(rows)))
    print("[OK] wrote {}".format(out))


if __name__ == "__main__":
    main()
