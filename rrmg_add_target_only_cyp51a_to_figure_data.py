#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rrmg_add_target_only_cyp51a_to_figure_data.py

Purpose:
  For each selected PlotID in figure_data/CYP51A_alignment/plot_manifest.tsv,
  map the target-only simulated Afum reads saved from Step 2 back to CYP51A,
  then copy/record these outputs beside existing figure_data.

This is designed for the comparison plot:
  - spiked metagenome reads -> CYP51A
  - Afum target-only reads  -> CYP51A

Typical HPC usage:
  cd /path/to/SpikeIN

  python3 rrmg_add_target_only_cyp51a_to_figure_data.py \
    --figure-dir figure_data/CYP51A_alignment \
    --gene ref/CYP51A_wt.fasta \
    --wd /path/to/SpikeIN \
    --threads 8

Output:
  figure_data/CYP51A_alignment/<PlotID>/target_only_reads_to_GENE.paf.gz
  figure_data/CYP51A_alignment/<PlotID>/target_only_reads_to_GENE.depth.tsv
  figure_data/CYP51A_alignment/<PlotID>/target_only_reads_to_GENE.summary.tsv
  figure_data/CYP51A_alignment/plot_manifest.with_target_only.tsv

Requirements:
  minimap2, samtools, gzip
"""

from __future__ import print_function

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Add target-only Afum reads-to-CYP51A PAF/depth outputs to RRMG figure_data."
    )
    p.add_argument("--figure-dir", default="figure_data/CYP51A_alignment",
                   help="Figure data directory containing plot_manifest.tsv.")
    p.add_argument("--manifest", default="",
                   help="Input manifest. Default: <figure-dir>/plot_manifest.tsv")
    p.add_argument("--out-manifest", default="",
                   help="Output manifest. Default: <figure-dir>/plot_manifest.with_target_only.tsv")
    p.add_argument("--gene", default="ref/CYP51A_wt.fasta",
                   help="Target gene FASTA, e.g. ref/CYP51A_wt.fasta")
    p.add_argument("--wd", default=".",
                   help="Project working directory.")
    p.add_argument("--root", default="temp/rrmg_permutation",
                   help="Permutation run root. Used if SourceRunDir is missing.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing outputs.")
    p.add_argument("--keep-bam", action="store_true",
                   help="Keep temporary target-only reads BAM and BAI.")
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--gzip", default="gzip")
    return p.parse_args()


def qrun(cmd, stdout_path=None, stderr_path=None, shell=False):
    if stdout_path:
        out_fh = open(stdout_path, "wb")
    else:
        out_fh = None

    if stderr_path:
        err_fh = open(stderr_path, "wb")
    else:
        err_fh = None

    try:
        subprocess.check_call(cmd, stdout=out_fh, stderr=err_fh, shell=shell)
    finally:
        if out_fh:
            out_fh.close()
        if err_fh:
            err_fh.close()


def resolve_path(p, wd):
    if not p:
        return ""
    p = str(p)
    if os.path.isabs(p):
        return p
    return str(Path(wd) / p)


def resolve_existing(p, wd):
    if not p:
        return ""

    candidates = []
    p = str(p)
    candidates.append(p)
    if not os.path.isabs(p):
        candidates.append(str(Path(wd) / p))

    # Imperial path variants
    candidates.append(p.replace("/rds/general/user/", "/rds/general/ephemeral/user/"))
    candidates.append(p.replace("/rds/general/ephemeral/user/", "/rds/general/user/"))

    for c in candidates:
        if c and os.path.exists(c):
            return c

    return candidates[0]


def safe(s):
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s).strip())
    return s.strip("_") or "NA"


def read_manifest(path):
    with open(path, newline="", errors="ignore") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_manifest(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def find_run_dir(row, wd, root):
    src = row.get("SourceRunDir", "")
    if src:
        p = resolve_existing(src, wd)
        if os.path.isdir(p):
            return Path(p)

    root_abs = Path(resolve_path(root, wd))
    rid = row.get("RunID", "")
    if rid:
        p = root_abs / safe(rid)
        if p.exists():
            return p

    seed = str(row.get("Seed", ""))
    ta = str(row.get("TargetAbundanceDesign", ""))
    pairs = str(row.get("TotalPairs", ""))

    if root_abs.exists():
        candidates = sorted(root_abs.glob("*Seed{}*".format(seed)))
        candidates = [
            c for c in candidates
            if ("ta{}".format(ta) in c.name and "pairs{}".format(pairs) in c.name)
        ]
        if candidates:
            return candidates[0]

    return None


def find_target_reads(run_dir, row):
    reads_dir = run_dir / "reads"
    if not reads_dir.exists():
        return None, None

    target_gid = safe(row.get("TargetGenomeID", "") or "Afum")

    candidates = sorted(reads_dir.glob("*target_1.fastq.gz"))
    preferred = [x for x in candidates if target_gid in x.name or "Afum" in x.name]
    if preferred:
        candidates = preferred

    if not candidates:
        return None, None

    r1 = candidates[0]
    r2 = Path(str(r1).replace("_1.fastq.gz", "_2.fastq.gz"))
    if not r2.exists():
        return None, None

    return r1, r2


def read_depth_summary(depth_file):
    depths = []
    with open(depth_file, errors="ignore") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                depths.append(int(float(parts[2])))
            except Exception:
                pass

    if not depths:
        return 0, 0, 0, 0

    covered = sum(1 for x in depths if x > 0)
    mean_depth = sum(depths) / float(len(depths))
    breadth = covered / float(len(depths)) * 100
    max_depth = max(depths)
    return covered, breadth, mean_depth, max_depth


def main():
    args = parse_args()

    wd = Path(args.wd).resolve()
    figure_dir = Path(resolve_path(args.figure_dir, wd))

    manifest = Path(args.manifest) if args.manifest else figure_dir / "plot_manifest.tsv"
    out_manifest = Path(args.out_manifest) if args.out_manifest else figure_dir / "plot_manifest.with_target_only.tsv"
    gene = Path(resolve_path(args.gene, wd))

    if not manifest.exists():
        raise SystemExit("[ERROR] manifest not found: {}".format(manifest))
    if not gene.exists():
        raise SystemExit("[ERROR] gene FASTA not found: {}".format(gene))

    rows = read_manifest(str(manifest))
    if not rows:
        raise SystemExit("[ERROR] empty manifest: {}".format(manifest))

    add_fields = [
        "TargetOnlyReadPAFGZ",
        "TargetOnlyReadDepthTSV",
        "TargetOnlyReadMappedAlignments",
        "TargetOnlyReadCoveredBp",
        "TargetOnlyReadBreadthPct",
        "TargetOnlyReadMeanDepth",
        "TargetOnlyReadMaxDepth",
        "TargetOnlyR1",
        "TargetOnlyR2",
    ]

    fieldnames = list(rows[0].keys())
    for f in add_fields:
        if f not in fieldnames:
            fieldnames.append(f)

    print("[INFO] manifest:", manifest)
    print("[INFO] rows:", len(rows))
    print("[INFO] gene:", gene)

    for i, row in enumerate(rows, 1):
        plot_id = row.get("PlotID", "") or "Plot{}".format(i)
        run_id = row.get("RunID", "")

        run_dir = find_run_dir(row, str(wd), args.root)
        if run_dir is None:
            print("[WARN] run_dir not found for PlotID={} RunID={}".format(plot_id, run_id))
            continue

        r1, r2 = find_target_reads(run_dir, row)
        if r1 is None:
            print("[WARN] target-only reads not found for PlotID={} RunID={}".format(plot_id, run_id))
            continue

        out_dir = figure_dir / plot_id
        out_dir.mkdir(parents=True, exist_ok=True)

        paf_gz = out_dir / "target_only_reads_to_GENE.paf.gz"
        depth_tsv = out_dir / "target_only_reads_to_GENE.depth.tsv"
        bam = out_dir / "target_only_reads_to_GENE.bam"
        summary_tsv = out_dir / "target_only_reads_to_GENE.summary.tsv"

        if args.force or not (paf_gz.exists() and paf_gz.stat().st_size > 0):
            print("[PAF] {}/{} {}".format(i, len(rows), plot_id))
            cmd = (
                "{minimap2} -t {threads} -x sr --secondary=no {gene} {r1} {r2} "
                "2> {log} | {gzip} -c > {out}"
            ).format(
                minimap2=args.minimap2,
                threads=args.threads,
                gene=str(gene),
                r1=str(r1),
                r2=str(r2),
                log=str(out_dir / "target_only_reads_to_GENE.paf.minimap2.log"),
                gzip=args.gzip,
                out=str(paf_gz)
            )
            qrun(cmd, shell=True)

        if args.force or not (depth_tsv.exists() and depth_tsv.stat().st_size > 0):
            print("[BAM/DEPTH] {}/{} {}".format(i, len(rows), plot_id))
            cmd = (
                "{minimap2} -t {threads} -ax sr --secondary=no {gene} {r1} {r2} "
                "2> {log} | "
                "{samtools} view -@ {threads} -b -F 4 - | "
                "{samtools} sort -@ {threads} -o {bam} -"
            ).format(
                minimap2=args.minimap2,
                threads=args.threads,
                gene=str(gene),
                r1=str(r1),
                r2=str(r2),
                log=str(out_dir / "target_only_reads_to_GENE.bam.minimap2.log"),
                samtools=args.samtools,
                bam=str(bam)
            )
            qrun(cmd, shell=True)

            qrun([args.samtools, "index", str(bam)])
            qrun([args.samtools, "depth", "-aa", str(bam)], stdout_path=str(depth_tsv))

        mapped_alignments = ""
        if bam.exists():
            try:
                mapped_alignments = subprocess.check_output([args.samtools, "view", "-c", str(bam)]).decode().strip()
            except Exception:
                mapped_alignments = ""

        covered_bp, breadth, mean_depth, max_depth = read_depth_summary(str(depth_tsv))

        with summary_tsv.open("w", newline="") as fh:
            fields = [
                "PlotID", "RunID", "SourceRunDir", "TargetOnlyR1", "TargetOnlyR2",
                "MappedAlignments", "CoveredBp", "BreadthPct", "MeanDepth", "MaxDepth",
                "PAFGZ", "DepthTSV"
            ]
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerow({
                "PlotID": plot_id,
                "RunID": run_id,
                "SourceRunDir": str(run_dir),
                "TargetOnlyR1": str(r1),
                "TargetOnlyR2": str(r2),
                "MappedAlignments": mapped_alignments,
                "CoveredBp": covered_bp,
                "BreadthPct": "{:.6f}".format(breadth),
                "MeanDepth": "{:.6f}".format(mean_depth),
                "MaxDepth": max_depth,
                "PAFGZ": str(paf_gz),
                "DepthTSV": str(depth_tsv),
            })

        row["TargetOnlyReadPAFGZ"] = "{}/{}".format(plot_id, paf_gz.name)
        row["TargetOnlyReadDepthTSV"] = "{}/{}".format(plot_id, depth_tsv.name)
        row["TargetOnlyReadMappedAlignments"] = mapped_alignments
        row["TargetOnlyReadCoveredBp"] = str(covered_bp)
        row["TargetOnlyReadBreadthPct"] = "{:.6f}".format(breadth)
        row["TargetOnlyReadMeanDepth"] = "{:.6f}".format(mean_depth)
        row["TargetOnlyReadMaxDepth"] = str(max_depth)
        row["TargetOnlyR1"] = str(r1)
        row["TargetOnlyR2"] = str(r2)

        if not args.keep_bam:
            for f in [bam, Path(str(bam) + ".bai")]:
                if f.exists():
                    f.unlink()

    write_manifest(out_manifest, rows, fieldnames)
    print("[OK] wrote:", out_manifest)


if __name__ == "__main__":
    main()
