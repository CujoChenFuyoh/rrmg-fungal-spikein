#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Evaluate one permutation's VAMB clusters against the known spike-in genome.

This evaluator is designed for simulated spike-in experiments where the target
reference genome is known.  It does not use CheckM/CheckM2/GUNC because their
quality estimates are intended for prokaryotic MAGs.  Instead, it reports:

  * target-aligned base pairs in each VAMB cluster;
  * target precision: target-aligned query bp / total cluster bp;
  * capture of the target-derived assembled contigs;
  * breadth of the target reference covered by each cluster.

The PAF input must be produced with target reference as minimap2 target and the
metagenomic contigs as query, for example:

  minimap2 -x asm5 target_genome.fasta contigs.fasta > contigs.to_target.paf
"""

from __future__ import print_function

import argparse
import csv
import os
import sys
from collections import defaultdict


SUMMARY_FIELDS = [
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
    "ContigsFA",
    "InputContigs",
    "InputBp",
    "VambAssignedContigs",
    "VambClusters",
    "EmittedBins",
    "TargetRef",
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
    "BAM",
    "ClustersTSV",
    "PerBinTSV",
    "VambOutDir",
]


PER_BIN_FIELDS = [
    "RankByTargetBp",
    "Bin",
    "BinFASTA",
    "Emitted",
    "BinBp",
    "NContigs",
    "TargetAlignedBp",
    "TargetPrecisionPct",
    "TargetAssemblyCapturePct",
    "TargetReferenceCoveredBp",
    "TargetReferenceBreadthPct",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Evaluate VAMB clusters against a known spike-in target genome."
    )
    p.add_argument("--clusters", required=True, help="VAMB cluster TSV.")
    p.add_argument("--contigs", required=True, help="Assembly contigs FASTA.")
    p.add_argument("--target-ref", required=True, help="Known spike-in genome FASTA.")
    p.add_argument("--target-paf", required=True, help="contigs-to-target minimap2 PAF.")
    p.add_argument("--bins-dir", default="", help="VAMB emitted FASTA bin directory.")
    p.add_argument("--out-bins", required=True, help="Output per-bin metrics TSV.")
    p.add_argument("--out-summary", required=True, help="Output one-row summary TSV.")
    p.add_argument("--run-index", default="")
    p.add_argument("--run-id", default="")
    p.add_argument("--sample", default="")
    p.add_argument("--seed", default="")
    p.add_argument("--target-abundance", default="")
    p.add_argument("--total-pairs", default="")
    p.add_argument("--target-expected-cov", default="")
    p.add_argument("--vamb-mode", default="")
    p.add_argument("--vamb-version", default="")
    p.add_argument("--vamb-seed", default="")
    p.add_argument("--min-contig", default="")
    p.add_argument("--min-bin-size", default="")
    p.add_argument("--bam", default="")
    p.add_argument("--vamb-outdir", default="")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-identity", type=float, default=0.95)
    p.add_argument("--min-aln-len", type=int, default=500)
    return p.parse_args()


def clean_cell(value):
    if value is None:
        return ""
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def fmt_pct(numerator, denominator):
    if not denominator:
        return "0.000000"
    return "{:.6f}".format(100.0 * float(numerator) / float(denominator))


def read_fasta_lengths(path):
    lengths = {}
    name = None
    length = 0
    with open(path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split(None, 1)[0]
                if not name:
                    raise ValueError("Empty FASTA identifier in {}".format(path))
                if name in lengths:
                    raise ValueError("Duplicate FASTA identifier: {}".format(name))
                length = 0
            else:
                if name is None:
                    raise ValueError("Sequence before FASTA header in {}".format(path))
                length += len("".join(line.split()))
    if name is not None:
        lengths[name] = length
    return lengths


def read_clusters(path):
    clusters = defaultdict(list)
    contig_to_cluster = {}
    with open(path, "rt", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row_number, row in enumerate(reader, 1):
            if not row or all(not x.strip() for x in row):
                continue
            if len(row) < 2:
                raise ValueError(
                    "Expected two columns in {}, row {}".format(path, row_number)
                )
            cluster = row[0].strip()
            contig = row[1].strip()
            if row_number == 1:
                left = cluster.lower().replace("_", "")
                right = contig.lower().replace("_", "")
                if "cluster" in left and ("contig" in right or "sequence" in right):
                    continue
            if not cluster or not contig:
                continue
            previous = contig_to_cluster.get(contig)
            if previous is not None and previous != cluster:
                raise ValueError(
                    "Contig {} occurs in multiple VAMB clusters".format(contig)
                )
            contig_to_cluster[contig] = cluster
            clusters[cluster].append(contig)
    return clusters, contig_to_cluster


def merge_intervals(intervals):
    if not intervals:
        return 0
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    start, end = ordered[0]
    for next_start, next_end in ordered[1:]:
        if next_start <= end:
            if next_end > end:
                end = next_end
        else:
            total += end - start
            start, end = next_start, next_end
    total += end - start
    return total


def read_filtered_paf(path, valid_queries, min_mapq, min_identity, min_aln_len):
    query_alignments = defaultdict(list)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return query_alignments
    with open(path, "rt") as fh:
        for row_number, raw in enumerate(fh, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError("Malformed PAF row {} in {}".format(row_number, path))
            query = fields[0]
            if query not in valid_queries:
                continue
            try:
                qstart = int(fields[2])
                qend = int(fields[3])
                target = fields[5]
                tstart = int(fields[7])
                tend = int(fields[8])
                matches = int(fields[9])
                block_len = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                raise ValueError(
                    "Non-numeric required PAF field at row {} in {}".format(
                        row_number, path
                    )
                )
            identity = float(matches) / block_len if block_len else 0.0
            if (
                mapq < min_mapq
                or block_len < min_aln_len
                or identity < min_identity
            ):
                continue
            query_alignments[query].append(
                (qstart, qend, target, tstart, tend)
            )
    return query_alignments


def target_metrics_for_contigs(contigs, query_alignments):
    query_bp = 0
    target_intervals = defaultdict(list)
    for contig in contigs:
        alignments = query_alignments.get(contig, [])
        query_bp += merge_intervals([(x[0], x[1]) for x in alignments])
        for _qstart, _qend, target, tstart, tend in alignments:
            target_intervals[target].append((tstart, tend))
    target_bp = sum(merge_intervals(x) for x in target_intervals.values())
    return query_bp, target_bp


def emitted_bin_index(bins_dir):
    index = {}
    if not bins_dir or not os.path.isdir(bins_dir):
        return index
    valid_extensions = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
    for name in sorted(os.listdir(bins_dir)):
        full = os.path.join(bins_dir, name)
        if not os.path.isfile(full):
            continue
        lower = name.lower()
        matched = None
        for extension in valid_extensions:
            if lower.endswith(extension):
                matched = extension
                break
        if matched is None:
            continue
        stem = name[: -len(matched)]
        index[stem] = os.path.abspath(full)
    return index


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
            writer.writerow({key: clean_cell(row.get(key, "")) for key in fields})


def main():
    args = parse_args()
    for path in (args.clusters, args.contigs, args.target_ref):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise SystemExit("[ERROR] Required non-empty file not found: {}".format(path))

    contig_lengths = read_fasta_lengths(args.contigs)
    target_lengths = read_fasta_lengths(args.target_ref)
    clusters, contig_to_cluster = read_clusters(args.clusters)
    if not clusters:
        raise SystemExit("[ERROR] No clusters parsed from {}".format(args.clusters))

    unknown_contigs = sorted(set(contig_to_cluster) - set(contig_lengths))
    if unknown_contigs:
        preview = ", ".join(unknown_contigs[:5])
        raise SystemExit(
            "[ERROR] VAMB cluster contigs absent from FASTA (first values: {})".format(
                preview
            )
        )

    query_alignments = read_filtered_paf(
        args.target_paf,
        set(contig_lengths),
        args.min_mapq,
        args.min_identity,
        args.min_aln_len,
    )
    all_target_query_bp, all_target_reference_bp = target_metrics_for_contigs(
        contig_lengths.keys(), query_alignments
    )
    target_genome_bp = sum(target_lengths.values())
    emitted = emitted_bin_index(args.bins_dir)

    metrics = []
    for cluster, contigs in clusters.items():
        bin_bp = sum(contig_lengths[x] for x in contigs)
        target_query_bp, target_reference_bp = target_metrics_for_contigs(
            contigs, query_alignments
        )
        bin_fasta = emitted.get(cluster, "")
        metrics.append(
            {
                "Bin": cluster,
                "BinFASTA": bin_fasta,
                "Emitted": "Y" if bin_fasta else "N",
                "BinBp": str(bin_bp),
                "NContigs": str(len(contigs)),
                "TargetAlignedBp": str(target_query_bp),
                "TargetPrecisionPct": fmt_pct(target_query_bp, bin_bp),
                "TargetAssemblyCapturePct": fmt_pct(
                    target_query_bp, all_target_query_bp
                ),
                "TargetReferenceCoveredBp": str(target_reference_bp),
                "TargetReferenceBreadthPct": fmt_pct(
                    target_reference_bp, target_genome_bp
                ),
                "_target_query_bp": target_query_bp,
                "_target_reference_bp": target_reference_bp,
                "_bin_bp": bin_bp,
            }
        )

    metrics.sort(
        key=lambda row: (
            -row["_target_query_bp"],
            -float(row["TargetPrecisionPct"]),
            -row["_target_reference_bp"],
            -row["_bin_bp"],
            row["Bin"],
        )
    )
    for rank, row in enumerate(metrics, 1):
        row["RankByTargetBp"] = str(rank)

    write_tsv(args.out_bins, PER_BIN_FIELDS, metrics)

    best = metrics[0] if metrics and metrics[0]["_target_query_bp"] > 0 else {}
    summary = {
        "RunIndex": args.run_index,
        "RunID": args.run_id,
        "Sample": args.sample,
        "Seed": args.seed,
        "TargetAbundance": args.target_abundance,
        "TotalPairs": args.total_pairs,
        "TargetExpectedCov": args.target_expected_cov,
        "Status": "ok",
        "VambMode": args.vamb_mode,
        "VambVersion": args.vamb_version,
        "VambSeed": args.vamb_seed,
        "MinContigBp": args.min_contig,
        "MinBinSizeBp": args.min_bin_size,
        "ContigsFA": os.path.abspath(args.contigs),
        "InputContigs": str(len(contig_lengths)),
        "InputBp": str(sum(contig_lengths.values())),
        "VambAssignedContigs": str(len(contig_to_cluster)),
        "VambClusters": str(len(clusters)),
        "EmittedBins": str(len(emitted)),
        "TargetRef": os.path.abspath(args.target_ref),
        "TargetGenomeBp": str(target_genome_bp),
        "TargetAlignedAssemblyBp": str(all_target_query_bp),
        "TargetAssemblyReferenceCoveredBp": str(all_target_reference_bp),
        "TargetAssemblyReferenceBreadthPct": fmt_pct(
            all_target_reference_bp, target_genome_bp
        ),
        "BestBinHasTargetAlignment": "Y" if best else "N",
        "BestBin": best.get("Bin", ""),
        "BestBinFASTA": best.get("BinFASTA", ""),
        "BestBinEmitted": best.get("Emitted", "N"),
        "BestBinBp": best.get("BinBp", "0"),
        "BestBinContigs": best.get("NContigs", "0"),
        "BestBinTargetAlignedBp": best.get("TargetAlignedBp", "0"),
        "BestBinTargetPrecisionPct": best.get("TargetPrecisionPct", "0.000000"),
        "BestBinTargetAssemblyCapturePct": best.get(
            "TargetAssemblyCapturePct", "0.000000"
        ),
        "BestBinTargetReferenceCoveredBp": best.get(
            "TargetReferenceCoveredBp", "0"
        ),
        "BestBinTargetReferenceBreadthPct": best.get(
            "TargetReferenceBreadthPct", "0.000000"
        ),
        "BAM": os.path.abspath(args.bam) if args.bam else "",
        "ClustersTSV": os.path.abspath(args.clusters),
        "PerBinTSV": os.path.abspath(args.out_bins),
        "VambOutDir": os.path.abspath(args.vamb_outdir)
        if args.vamb_outdir
        else "",
    }
    write_tsv(args.out_summary, SUMMARY_FIELDS, [summary])
    print("[OK] Per-bin VAMB metrics: {}".format(args.out_bins))
    print("[OK] Per-run VAMB summary: {}".format(args.out_summary))


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        pass
    except Exception as exc:
        print("[ERROR] {}".format(exc), file=sys.stderr)
        raise
