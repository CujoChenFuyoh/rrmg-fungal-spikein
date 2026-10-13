#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Collect one-permutation VAMB summaries and aggregate recovery by condition.

The permutation plan is used as the authoritative run list, so failed or
unsubmitted array rows appear as Status=missing instead of silently vanishing.
"""

from __future__ import print_function

import argparse
import csv
import os
import re
import statistics
from collections import defaultdict


PREFERRED_FIELDS = [
    "RunIndex",
    "RunID",
    "Sample",
    "Seed",
    "TargetAbundance",
    "TotalPairs",
    "TargetExpectedCov",
    "Status",
    "VambMode",
    "VambVersion",
    "VambSeed",
    "MinContigBp",
    "MinBinSizeBp",
    "InputContigs",
    "InputBp",
    "VambAssignedContigs",
    "VambClusters",
    "EmittedBins",
    "TargetGenomeBp",
    "TargetAlignedAssemblyBp",
    "TargetAssemblyReferenceCoveredBp",
    "TargetAssemblyReferenceBreadthPct",
    "BestBinHasTargetAlignment",
    "BestBin",
    "BestBinFASTA",
    "BestBinEmitted",
    "BestBinBp",
    "BestBinContigs",
    "BestBinTargetAlignedBp",
    "BestBinTargetPrecisionPct",
    "BestBinTargetAssemblyCapturePct",
    "BestBinTargetReferenceCoveredBp",
    "BestBinTargetReferenceBreadthPct",
    "ContigsFA",
    "TargetRef",
    "BAM",
    "ClustersTSV",
    "PerBinTSV",
    "VambOutDir",
    "SummaryTSV",
]


AGGREGATE_FIELDS = [
    "TargetAbundance",
    "TotalPairs",
    "Runs",
    "CompletedRuns",
    "CompletionRatePct",
    "RunsWithTargetAlignment",
    "TargetAlignmentRatePct",
    "RunsWithEmittedTargetBin",
    "EmittedTargetBinRatePct",
    "MeanTargetAssemblyReferenceBreadthPct",
    "MedianTargetAssemblyReferenceBreadthPct",
    "MeanBestBinTargetPrecisionPct",
    "MedianBestBinTargetPrecisionPct",
    "MeanBestBinTargetAssemblyCapturePct",
    "MedianBestBinTargetAssemblyCapturePct",
    "MeanBestBinTargetReferenceBreadthPct",
    "MedianBestBinTargetReferenceBreadthPct",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Collect VAMB per-run summaries for the RRMG permutation workflow."
    )
    p.add_argument("--plan", required=True, help="rrmg_permutation_plan.tsv")
    p.add_argument(
        "--root",
        default="temp/rrmg_permutation",
        help="Root containing one directory per run.",
    )
    p.add_argument(
        "--out",
        default="result/rrmg_vamb_per_run.tsv",
        help="Detailed output TSV.",
    )
    p.add_argument(
        "--wd",
        default=".",
        help="Working directory used to resolve relative paths.",
    )
    p.add_argument(
        "--no-aggregate",
        action="store_true",
        help="Do not create .by_abundance_depth.tsv.",
    )
    return p.parse_args()


def safe_tag(value):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return value.strip("_") or "NA"


def resolve(path, wd):
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(wd, path))


def read_tsv(path):
    with open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def pct(numerator, denominator):
    if not denominator:
        return "0.000000"
    return "{:.6f}".format(100.0 * float(numerator) / float(denominator))


def mean_or_blank(values):
    values = [x for x in values if x is not None]
    if not values:
        return ""
    return "{:.6f}".format(sum(values) / float(len(values)))


def median_or_blank(values):
    values = [x for x in values if x is not None]
    if not values:
        return ""
    return "{:.6f}".format(statistics.median(values))


def output_aggregate_path(out_path):
    root, extension = os.path.splitext(out_path)
    if extension.lower() == ".tsv":
        return root + ".by_abundance_depth.tsv"
    return out_path + ".by_abundance_depth.tsv"


def write_tsv(path, fields, rows):
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)
    with open(path, "wt", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def main():
    args = parse_args()
    wd = os.path.abspath(args.wd)
    plan_path = resolve(args.plan, wd)
    root = resolve(args.root, wd)
    out = resolve(args.out, wd)

    plan_rows = read_tsv(plan_path)
    detailed = []
    extra_fields = []
    seen_extra = set()

    for plan_row in plan_rows:
        run_index = plan_row.get("run_index", "")
        run_id = plan_row.get("run_id") or plan_row.get("sample") or "row{}".format(
            run_index
        )
        sample = plan_row.get("sample") or run_id
        seed = plan_row.get("seed", "")
        run_dir = os.path.join(root, safe_tag(run_id))
        summary_path = os.path.join(run_dir, "binning", "vamb_summary.tsv")

        summary_rows = read_tsv(summary_path) if os.path.isfile(summary_path) else []
        if summary_rows:
            row = dict(summary_rows[0])
            row["SummaryTSV"] = os.path.abspath(summary_path)
            for key in row:
                if key not in PREFERRED_FIELDS and key not in seen_extra:
                    seen_extra.add(key)
                    extra_fields.append(key)
        else:
            row = {
                "RunIndex": run_index,
                "RunID": run_id,
                "Sample": sample,
                "Seed": seed,
                "TargetAbundance": plan_row.get("target_abundance", ""),
                "TotalPairs": plan_row.get("total_pairs", ""),
                "TargetExpectedCov": plan_row.get("target_expected_cov", ""),
                "Status": "missing",
                "SummaryTSV": os.path.abspath(summary_path),
            }
        detailed.append(row)

    fields = PREFERRED_FIELDS + extra_fields
    write_tsv(out, fields, detailed)
    print("[OK] Detailed VAMB table: {}".format(out))

    if args.no_aggregate:
        return

    groups = defaultdict(list)
    for row in detailed:
        groups[(row.get("TargetAbundance", ""), row.get("TotalPairs", ""))].append(
            row
        )

    aggregate_rows = []
    for key in sorted(
        groups,
        key=lambda x: (
            as_float(x[0]) if as_float(x[0]) is not None else float("inf"),
            as_float(x[1]) if as_float(x[1]) is not None else float("inf"),
            x,
        ),
    ):
        rows = groups[key]
        completed = [row for row in rows if row.get("Status") == "ok"]
        target_aligned = [
            row
            for row in completed
            if row.get("BestBinHasTargetAlignment", "").upper() == "Y"
        ]
        emitted = [
            row
            for row in target_aligned
            if row.get("BestBinEmitted", "").upper() == "Y"
        ]

        assembly_breadth = [
            as_float(row.get("TargetAssemblyReferenceBreadthPct"))
            for row in completed
        ]
        best_precision = [
            as_float(row.get("BestBinTargetPrecisionPct")) for row in completed
        ]
        best_capture = [
            as_float(row.get("BestBinTargetAssemblyCapturePct"))
            for row in completed
        ]
        best_breadth = [
            as_float(row.get("BestBinTargetReferenceBreadthPct"))
            for row in completed
        ]

        aggregate_rows.append(
            {
                "TargetAbundance": key[0],
                "TotalPairs": key[1],
                "Runs": str(len(rows)),
                "CompletedRuns": str(len(completed)),
                "CompletionRatePct": pct(len(completed), len(rows)),
                "RunsWithTargetAlignment": str(len(target_aligned)),
                "TargetAlignmentRatePct": pct(len(target_aligned), len(completed)),
                "RunsWithEmittedTargetBin": str(len(emitted)),
                "EmittedTargetBinRatePct": pct(len(emitted), len(completed)),
                "MeanTargetAssemblyReferenceBreadthPct": mean_or_blank(
                    assembly_breadth
                ),
                "MedianTargetAssemblyReferenceBreadthPct": median_or_blank(
                    assembly_breadth
                ),
                "MeanBestBinTargetPrecisionPct": mean_or_blank(best_precision),
                "MedianBestBinTargetPrecisionPct": median_or_blank(best_precision),
                "MeanBestBinTargetAssemblyCapturePct": mean_or_blank(best_capture),
                "MedianBestBinTargetAssemblyCapturePct": median_or_blank(best_capture),
                "MeanBestBinTargetReferenceBreadthPct": mean_or_blank(best_breadth),
                "MedianBestBinTargetReferenceBreadthPct": median_or_blank(
                    best_breadth
                ),
            }
        )

    aggregate_path = output_aggregate_path(out)
    write_tsv(aggregate_path, AGGREGATE_FIELDS, aggregate_rows)
    print("[OK] Aggregated VAMB table: {}".format(aggregate_path))


if __name__ == "__main__":
    main()
