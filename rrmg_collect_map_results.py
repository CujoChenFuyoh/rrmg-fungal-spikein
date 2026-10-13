#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rrmg_collect_map_results.py

Merge RRMG split-pipeline summary, plan, per-run map/target_call_summary.tsv,
and per-run design/merge_design.tsv into one target-gene detection table.

This script is intentionally dependency-free and should run with Python 3.6+.

Typical usage:
  python3 rrmg_collect_map_results.py \
    --summary result/rrmg_split_pipeline_summary.tsv \
    --plan input/rrmg/rrmg_permutation_plan.tsv \
    --root temp/rrmg_permutation \
    --out result/rrmg_target_gene_map_merged.tsv

Outputs:
  1. Detailed per-run table:
     result/rrmg_target_gene_map_merged.tsv

  2. Aggregated hit-rate table:
     result/rrmg_target_gene_map_merged.by_abundance_depth.tsv

Interpretation note:
  In the current split workflow, HasHit/identity/breadth describe detection of
  the target gene. The Step 3 VCF is only a placeholder, so VarN must not be
  interpreted as a validated fungAMR mutation count.
"""

from __future__ import print_function

import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge RRMG map target-gene hit results with abundance/depth metadata."
    )
    p.add_argument(
        "--summary",
        required=True,
        help="Main split-pipeline summary TSV, e.g. result/rrmg_split_pipeline_summary.tsv",
    )
    p.add_argument(
        "--plan",
        default="",
        help="Permutation plan TSV, e.g. input/rrmg/rrmg_permutation_plan.tsv",
    )
    p.add_argument(
        "--root",
        default="temp/rrmg_permutation",
        help="Root directory containing run dirs. Default: temp/rrmg_permutation",
    )
    p.add_argument(
        "--out",
        default="result/rrmg_target_gene_map_merged.tsv",
        help="Output detailed merged TSV.",
    )
    p.add_argument(
        "--readlen",
        type=int,
        default=150,
        help="Default read length if not found in design/plan. Default: 150",
    )
    p.add_argument(
        "--wd",
        default=".",
        help="Working directory used to resolve relative paths. Default: current directory.",
    )
    p.add_argument(
        "--no-aggregate",
        action="store_true",
        help="Do not write the by-abundance/depth hit-rate summary.",
    )
    return p.parse_args()


def safe_str(x):
    if x is None:
        return ""
    return str(x)


def is_file(path):
    return bool(path) and os.path.isfile(path)


def resolve_path(path, wd):
    if not path:
        return ""
    path = str(path)
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(wd, path))


def try_resolve_existing(path, wd):
    """Resolve a path and return it if it exists, otherwise return original resolved path."""
    if not path:
        return ""
    candidates = []

    if os.path.isabs(path):
        candidates.append(path)
    else:
        candidates.append(os.path.join(wd, path))
        candidates.append(path)

    # Handle common Imperial RDS path variants seen in logs.
    extra = []
    for c in candidates:
        extra.append(c.replace("/rds/general/user/", "/rds/general/ephemeral/user/"))
        extra.append(c.replace("/rds/general/ephemeral/user/", "/rds/general/user/"))
    candidates.extend(extra)

    for c in candidates:
        if os.path.exists(c):
            return os.path.abspath(c)

    return os.path.abspath(candidates[0])


def read_tsv(path):
    path = try_resolve_existing(path, os.getcwd())
    if not os.path.isfile(path):
        raise SystemExit("[ERROR] TSV not found: {}".format(path))
    with open(path, newline="", errors="ignore") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_key_value_tsv(path):
    """Read key/value TSV like design/merge_design.tsv."""
    d = {}
    if not path or not os.path.isfile(path):
        return d
    with open(path, errors="ignore") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if len(header) >= 2 and header[0].lower() == "key":
            for line in fh:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) == 2:
                    d[parts[0]] = parts[1]
        else:
            # Fallback: parse all two-column lines.
            fh.seek(0)
            for line in fh:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) == 2:
                    d[parts[0]] = parts[1]
    return d


def open_text(path):
    if not path:
        raise IOError("empty path")
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", errors="ignore")


def fasta_len(path):
    path = try_resolve_existing(path, os.getcwd())
    if not os.path.isfile(path):
        return ""
    n = 0
    try:
        with open_text(path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    n += len(line.strip())
    except Exception:
        return ""
    return str(n)


def to_float(x, default=None):
    if x is None or str(x) == "":
        return default
    try:
        return float(x)
    except Exception:
        return default


def to_int(x, default=None):
    if x is None or str(x) == "":
        return default
    try:
        return int(round(float(x)))
    except Exception:
        return default


def format_float(x, nd=10):
    if x is None:
        return ""
    try:
        return ("{0:." + str(nd) + "f}").format(float(x))
    except Exception:
        return safe_str(x)


def derive_run_id_from_dir(path):
    if not path:
        return ""
    path = path.rstrip("/")
    return os.path.basename(path)


def build_map_index(root, wd):
    root = try_resolve_existing(root, wd)
    idx = {}
    if not os.path.isdir(root):
        return idx

    for dirpath, dirnames, filenames in os.walk(root):
        if "target_call_summary.tsv" in filenames and os.path.basename(dirpath) == "map":
            f = os.path.join(dirpath, "target_call_summary.tsv")
            run_dir = os.path.dirname(dirpath)
            run_id = os.path.basename(run_dir)
            idx[run_id] = f
    return idx


def find_existing_map(summary_row, map_index, root, wd):
    # 1. Direct MapSummary path
    ms = summary_row.get("MapSummary", "")
    ms2 = try_resolve_existing(ms, wd) if ms else ""
    if is_file(ms2):
        return ms2

    # 2. RunDir/map/target_call_summary.tsv
    rd = summary_row.get("RunDir", "")
    rd2 = try_resolve_existing(rd, wd) if rd else ""
    cand = os.path.join(rd2, "map", "target_call_summary.tsv") if rd2 else ""
    if is_file(cand):
        return cand

    # 3. Use run dir basename in root index
    run_id = derive_run_id_from_dir(rd)
    if run_id and run_id in map_index:
        return map_index[run_id]

    # 4. Try reconstruct from known naming pattern.
    seed = summary_row.get("Seed", "")
    ta = summary_row.get("TargetAbundance", "")
    tp = summary_row.get("TotalPairs", "")
    if seed and ta and tp:
        maybe = "AirMG_Afum__ta{}__pairs{}__Seed{}".format(ta, tp, seed)
        if maybe in map_index:
            return map_index[maybe]

    return ""


def find_design(summary_row, map_path, root, wd):
    candidates = []

    rd = summary_row.get("RunDir", "")
    if rd:
        rd2 = try_resolve_existing(rd, wd)
        candidates.append(os.path.join(rd2, "design", "merge_design.tsv"))

    if map_path:
        run_dir = os.path.dirname(os.path.dirname(map_path))
        candidates.append(os.path.join(run_dir, "design", "merge_design.tsv"))

    for c in candidates:
        if os.path.isfile(c):
            return c
    return ""


def build_plan_index(plan_path, wd):
    by_run_id = {}
    by_seed_sample = {}
    if not plan_path:
        return by_run_id, by_seed_sample

    plan_path = try_resolve_existing(plan_path, wd)
    if not os.path.isfile(plan_path):
        print("[WARN] plan not found, continue without plan: {}".format(plan_path), file=sys.stderr)
        return by_run_id, by_seed_sample

    with open(plan_path, newline="", errors="ignore") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rid = row.get("run_id", "")
            if rid:
                by_run_id[rid] = row
            key = (str(row.get("seed", "")), str(row.get("sample", "")))
            by_seed_sample[key] = row

    return by_run_id, by_seed_sample


def read_first_map_row(path):
    if not path or not os.path.isfile(path):
        return {}
    with open(path, newline="", errors="ignore") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return {}
    # Current pipeline has one target gene row. If multiple rows exist, use first;
    # users can modify this to emit one row per target if needed.
    return rows[0]


def main():
    args = parse_args()
    wd = os.path.abspath(args.wd)
    summary_path = try_resolve_existing(args.summary, wd)
    root = try_resolve_existing(args.root, wd)

    summary = read_tsv(summary_path)
    plan_by_run, plan_by_seed_sample = build_plan_index(args.plan, wd)
    map_index = build_map_index(root, wd)

    if not summary:
        raise SystemExit("[ERROR] Empty summary: {}".format(summary_path))

    print("[INFO] summary rows: {}".format(len(summary)), file=sys.stderr)
    print("[INFO] map files indexed under root: {}".format(len(map_index)), file=sys.stderr)

    out_rows = []

    for srow in summary:
        map_path = find_existing_map(srow, map_index, root, wd)
        mrow = read_first_map_row(map_path)
        design_path = find_design(srow, map_path, root, wd)
        design = read_key_value_tsv(design_path)

        run_dir = srow.get("RunDir", "")
        run_id = design.get("run_id", "") or derive_run_id_from_dir(run_dir)
        if not run_id and map_path:
            run_id = os.path.basename(os.path.dirname(os.path.dirname(map_path)))

        seed = srow.get("Seed", "") or design.get("target_seed", "")
        sample = srow.get("Sample", "") or design.get("sample", "")

        prow = {}
        if run_id in plan_by_run:
            prow = plan_by_run[run_id]
        elif (str(seed), str(sample)) in plan_by_seed_sample:
            prow = plan_by_seed_sample[(str(seed), str(sample))]

        # Basic design fields
        target_abund = (
            design.get("target_abundance", "")
            or srow.get("TargetAbundance", "")
            or prow.get("target_abundance", "")
        )
        total_pairs = (
            design.get("total_pairs", "")
            or srow.get("TotalPairs", "")
            or prow.get("total_pairs", "")
        )
        readlen = design.get("readlen", "") or prow.get("readlen", "") or str(args.readlen)

        target_pairs = design.get("target_pairs", "")
        if not target_pairs:
            tp = to_float(total_pairs)
            ta = to_float(target_abund)
            if tp is not None and ta is not None:
                target_pairs = str(int(round(tp * ta)))

        target_genome_bp = design.get("target_genome_bp", "")
        target_fasta = design.get("target_fasta", "") or prow.get("target_fasta", "")
        if (not target_genome_bp) and target_fasta:
            target_genome_bp = fasta_len(target_fasta)

        target_cov = design.get("target_cov", "")
        if not target_cov:
            pairs = to_float(target_pairs)
            rl = to_float(readlen)
            gbp = to_float(target_genome_bp)
            if pairs is not None and rl is not None and gbp and gbp > 0:
                target_cov = format_float(pairs * rl * 2.0 / gbp, 10)

        target_expected_bases = ""
        pairs = to_float(target_pairs)
        rl = to_float(readlen)
        if pairs is not None and rl is not None:
            target_expected_bases = str(int(round(pairs * rl * 2.0)))

        # Target gene length
        target_fa = mrow.get("TargetFA", "")
        target_gene_bp = fasta_len(target_fa) if target_fa else ""

        covered_bp = mrow.get("TargetCoveredBp", "")
        target_breadth = mrow.get("TargetBreadthPct", "")

        uncovered_bp = ""
        if target_gene_bp and covered_bp:
            try:
                uncovered_bp = str(max(0, int(float(target_gene_bp)) - int(float(covered_bp))))
            except Exception:
                uncovered_bp = ""

        has_hit = mrow.get("HasHit", "")
        if has_hit == "":
            has_hit = "NA" if not map_path else "N"

        out = {
            "RunID": run_id,
            "Seed": seed,
            "Sample": sample,
            "TargetName": design.get("target_name", "") or prow.get("target_name", "") or srow.get("KrakenSpecies", ""),
            "TargetTaxID": design.get("target_taxid", "") or prow.get("target_taxid", "") or srow.get("SpeciesTaxID", ""),
            "TargetGenomeID": design.get("target_genome_id", "") or prow.get("target_genome_id", ""),
            "TargetAbundanceDesign": target_abund,
            "TotalPairs": total_pairs,
            "ReadLen": readlen,
            "TargetPairsExpected": target_pairs,
            "TargetExpectedBases": target_expected_bases,
            "TargetGenomeBp": target_genome_bp,
            "TargetGenomeCoverageExpected": target_cov,

            "KrakenSpecies": srow.get("KrakenSpecies", ""),
            "KrakenSpeciesPct": srow.get("SpeciesPct", ""),
            "KrakenSpeciesCladeReads": srow.get("SpeciesCladeReads", ""),
            "KrakenSpeciesTaxonReads": srow.get("SpeciesTaxonReads", ""),
            "KrakenSpeciesTaxID": srow.get("SpeciesTaxID", ""),

            "HasGeneHit": has_hit,
            "TargetType": mrow.get("TargetType", ""),
            "TargetFA": target_fa,
            "TargetGeneBp": target_gene_bp,
            "TargetBreadthPct": target_breadth,
            "TargetCoveredBp": covered_bp,
            "TargetUncoveredBp": uncovered_bp,
            "BestAlnLen": mrow.get("BestAlnLen", ""),
            "BestIdentity": mrow.get("BestIdentity", ""),
            "VarN": mrow.get("VarN", ""),
            "ContigsFA": mrow.get("ContigsFA", ""),
            "VCF": mrow.get("VCF", ""),

            "MapSummary": map_path,
            "MergeDesign": design_path,
            "RunDir": run_dir,
        }
        out_rows.append(out)

    fields = [
        "RunID", "Seed", "Sample",
        "TargetName", "TargetTaxID", "TargetGenomeID",
        "TargetAbundanceDesign", "TotalPairs", "ReadLen",
        "TargetPairsExpected", "TargetExpectedBases",
        "TargetGenomeBp", "TargetGenomeCoverageExpected",
        "KrakenSpecies", "KrakenSpeciesPct",
        "KrakenSpeciesCladeReads", "KrakenSpeciesTaxonReads", "KrakenSpeciesTaxID",
        "HasGeneHit", "TargetType", "TargetFA", "TargetGeneBp",
        "TargetBreadthPct", "TargetCoveredBp", "TargetUncoveredBp",
        "BestAlnLen", "BestIdentity", "VarN",
        "ContigsFA", "VCF", "MapSummary", "MergeDesign", "RunDir",
    ]

    out_path = resolve_path(args.out, wd)
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    print("[OK] wrote detail table: {}".format(out_path), file=sys.stderr)

    if not args.no_aggregate:
        agg = defaultdict(list)
        for r in out_rows:
            key = (r.get("TargetAbundanceDesign", ""), r.get("TotalPairs", ""))
            agg[key].append(r)

        agg_fields = [
            "TargetAbundanceDesign", "TotalPairs", "N",
            "HitN", "HitRate",
            "MedianTargetBreadthPct", "MeanTargetBreadthPct",
            "MedianKrakenSpeciesPct", "MeanKrakenSpeciesPct",
            "MedianTargetGenomeCoverageExpected",
        ]

        def median(vals):
            vals = sorted([v for v in vals if v is not None])
            if not vals:
                return ""
            n = len(vals)
            if n % 2:
                return vals[n // 2]
            return (vals[n // 2 - 1] + vals[n // 2]) / 2.0

        def mean(vals):
            vals = [v for v in vals if v is not None]
            if not vals:
                return ""
            return sum(vals) / float(len(vals))

        agg_path = out_path.rsplit(".", 1)[0] + ".by_abundance_depth.tsv"
        with open(agg_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=agg_fields, delimiter="\t", lineterminator="\n")
            w.writeheader()

            for key in sorted(agg, key=lambda x: (to_float(x[0], 0), to_float(x[1], 0))):
                rows = agg[key]
                n = len(rows)
                hit_n = sum(1 for r in rows if str(r.get("HasGeneHit", "")).upper() == "Y")
                breadth = [to_float(r.get("TargetBreadthPct")) for r in rows]
                kraken_pct = [to_float(r.get("KrakenSpeciesPct")) for r in rows]
                cov = [to_float(r.get("TargetGenomeCoverageExpected")) for r in rows]
                w.writerow({
                    "TargetAbundanceDesign": key[0],
                    "TotalPairs": key[1],
                    "N": n,
                    "HitN": hit_n,
                    "HitRate": format_float(hit_n / float(n) if n else 0, 4),
                    "MedianTargetBreadthPct": format_float(median(breadth), 6),
                    "MeanTargetBreadthPct": format_float(mean(breadth), 6),
                    "MedianKrakenSpeciesPct": format_float(median(kraken_pct), 6),
                    "MeanKrakenSpeciesPct": format_float(mean(kraken_pct), 6),
                    "MedianTargetGenomeCoverageExpected": format_float(median(cov), 10),
                })

        print("[OK] wrote aggregate table: {}".format(agg_path), file=sys.stderr)


if __name__ == "__main__":
    main()
