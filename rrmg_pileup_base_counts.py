#!/usr/bin/env python3
"""Convert one CYP51A BAM into an auditable per-position base-count table.

The script calls ``samtools mpileup -aa`` and parses every mapped base after
the requested mapping/base-quality filters.  It records strand-specific A/C/G/T
counts, SNV non-reference support and exact insertion/deletion event strings.
It does not make a mutation call; downstream code can apply predeclared DP/AD/
allele-fraction rules and join the result to the validated FungAMR site table.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile


BASES = "ACGTN"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--bam", type=Path)
    source.add_argument("--empty", action="store_true", help="write explicit zero counts for every reference base")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--label", required=True, help="metagenomic or synthetic_target")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--max-depth", type=int, default=0, help="0 removes the mpileup depth cap")
    return parser.parse_args()


def reference_inventory(path: Path) -> tuple[dict[str, str], str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                tokens = line[1:].split()
                name = tokens[0] if tokens else ""
                if not name or name in records:
                    raise ValueError("Reference FASTA has an empty or duplicate sequence name")
                records[name] = []
            elif name is None:
                raise ValueError("Reference FASTA sequence occurs before its header")
            else:
                seq = "".join(line.split()).upper()
                if any(base not in "ACGTRYSWKMBDHVN" for base in seq):
                    raise ValueError("Reference FASTA contains an invalid base")
                records[name].append(seq)
    joined = {key: "".join(parts) for key, parts in records.items()}
    if not joined or any(not seq for seq in joined.values()):
        raise ValueError("Reference FASTA is empty")
    if len(joined) != 1:
        raise ValueError("CYP51A base counts require exactly one reference sequence")
    sequence = next(iter(joined.values()))
    return joined, hashlib.sha256(sequence.encode()).hexdigest()


def add_base(counts: Counter[str], ref: str, symbol: str) -> None:
    if symbol == ".":
        base, strand = ref, "F"
    elif symbol == ",":
        base, strand = ref, "R"
    else:
        base = symbol.upper()
        strand = "F" if symbol.isupper() else "R"
    if base not in BASES:
        base = "N"
    counts[base] += 1
    counts[f"{base}_{strand}"] += 1


def parse_read_bases(ref: str, text: str, declared_depth: int) -> dict[str, object]:
    """Parse the read-bases column documented for samtools mpileup."""
    counts: Counter[str] = Counter()
    insertions: Counter[str] = Counter()
    deletions: Counter[str] = Counter()
    i = 0
    if declared_depth == 0:
        text = ""
    while i < len(text):
        symbol = text[i]
        if symbol == "^":
            if i + 1 >= len(text):
                raise ValueError("Truncated mpileup read-start marker")
            i += 2  # skip '^' and encoded mapping-quality character
        elif symbol == "$":
            i += 1
        elif symbol in ".," or symbol in "ACGTNacgtn":
            add_base(counts, ref, symbol)
            i += 1
        elif symbol in "+-":
            sign = symbol
            i += 1
            start = i
            while i < len(text) and text[i].isdigit():
                i += 1
            if start == i:
                raise ValueError("Mpileup indel has no length")
            length = int(text[start:i])
            allele = text[i : i + length].upper()
            if len(allele) != length:
                raise ValueError("Truncated mpileup indel allele")
            if sign == "+":
                insertions[allele] += 1
            else:
                deletions[allele] += 1
            i += length
        elif symbol in "*#":
            counts["DeletionPlaceholder"] += 1
            i += 1
        elif symbol in "<>":
            counts["ReferenceSkip"] += 1
            i += 1
        else:
            raise ValueError(f"Unsupported mpileup read-base symbol: {symbol!r}")

    ref = ref if ref in BASES else "N"
    ref_count = counts[ref]
    base_depth = sum(counts[base] for base in BASES)
    if base_depth + counts["DeletionPlaceholder"] + counts["ReferenceSkip"] != declared_depth:
        raise ValueError("Parsed read-base count differs from declared mpileup depth")
    nonref = sum(counts[base] for base in "ACGT" if base != ref)
    alternatives = [(counts[base], base) for base in "ACGT" if base != ref]
    best_count, best_base = max(alternatives, default=(0, "."))
    if best_count == 0:
        best_base = "."
    return {
        **{base: counts[base] for base in BASES},
        **{f"{base}_{strand}": counts[f"{base}_{strand}"] for base in BASES for strand in ("F", "R")},
        "BaseCountDepth": base_depth,
        "RefCount": ref_count,
        "SNVNonRefCount": nonref,
        "SNVNonRefFraction": (nonref / base_depth) if base_depth else 0.0,
        "MajorNonRefBase": best_base,
        "MajorNonRefCount": best_count,
        "InsertionEvents": sum(insertions.values()),
        "InsertionAlleles": ";".join(f"{key}:{value}" for key, value in sorted(insertions.items())),
        "DeletionStartEvents": sum(deletions.values()),
        "DeletionAlleles": ";".join(f"{key}:{value}" for key, value in sorted(deletions.items())),
        "DeletionPlaceholders": counts["DeletionPlaceholder"],
        "ReferenceSkips": counts["ReferenceSkip"],
        "ObservationStatus": "NO_COVERAGE" if declared_depth == 0 else
            "NO_ACGT_BASES" if not sum(counts[b] for b in "ACGT") else
            "NONREF_OBSERVED" if nonref else "REF_ONLY",
    }


FIELDS = [
    "Source", "Reference", "Position", "RefBase", "PileupDepth", "BaseCountDepth",
    "A", "C", "G", "T", "N",
    "A_F", "A_R", "C_F", "C_R", "G_F", "G_R", "T_F", "T_R", "N_F", "N_R",
    "RefCount", "SNVNonRefCount", "SNVNonRefFraction", "MajorNonRefBase", "MajorNonRefCount",
    "InsertionEvents", "InsertionAlleles", "DeletionStartEvents", "DeletionAlleles",
    "DeletionPlaceholders", "ReferenceSkips", "ObservationStatus", "MinMapQ", "MinBaseQ",
    "ReferenceSHA256", "BAM",
]


def output_handle(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("w", encoding="utf-8", newline="")


def main() -> int:
    args = parse_args()
    if args.min_mapq < 0 or args.min_baseq < 0 or args.max_depth < 0:
        raise SystemExit("[ERROR] Quality thresholds and max depth must be non-negative")
    if not args.reference.is_file() or not args.reference.stat().st_size:
        raise SystemExit(f"[ERROR] Missing/empty reference: {args.reference}")
    if args.bam is not None and (not args.bam.is_file() or not args.bam.stat().st_size):
        raise SystemExit(f"[ERROR] Missing/empty BAM: {args.bam}")
    references, ref_sha = reference_inventory(args.reference)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if str(args.output).endswith(".gz"):
        tmp = args.output.with_name(args.output.name[:-3] + f".tmp.{os.getpid()}.gz")
    else:
        tmp = args.output.with_name(args.output.name + f".tmp.{os.getpid()}")
    if args.empty:
        try:
            with output_handle(tmp) as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for chrom, sequence in references.items():
                    for position, ref in enumerate(sequence, 1):
                        evidence = parse_read_bases(ref, "", 0)
                        writer.writerow({
                            "Source": args.label, "Reference": chrom, "Position": position,
                            "RefBase": ref, "PileupDepth": 0, **evidence,
                            "SNVNonRefFraction": "0.00000000", "MinMapQ": args.min_mapq,
                            "MinBaseQ": args.min_baseq, "ReferenceSHA256": ref_sha, "BAM": "",
                        })
            tmp.replace(args.output)
            print(f"[DONE] {args.output} source={args.label} explicit_zero=1", flush=True)
            return 0
        except Exception:
            tmp.unlink(missing_ok=True)
            raise

    assert args.bam is not None
    command = [
        args.samtools, "mpileup", "-aa", "-A", "-B", "-x", "--ff", "0xF04", "-d", str(args.max_depth),
        "-q", str(args.min_mapq), "-Q", str(args.min_baseq),
        "-f", str(args.reference), str(args.bam),
    ]
    errors = tempfile.TemporaryFile(mode="w+t", encoding="utf-8")
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=errors,
        text=True, encoding="utf-8", errors="replace",
    )
    assert process.stdout is not None
    observed: Counter[str] = Counter()
    try:
        with output_handle(tmp) as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for raw in process.stdout:
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 6:
                    raise ValueError("samtools mpileup returned fewer than six columns")
                chrom, pos_text, ref, depth_text, read_bases = fields[:5]
                if chrom not in references:
                    raise ValueError(f"Mpileup sequence is absent from reference: {chrom}")
                position = int(pos_text)
                depth = int(depth_text)
                if not 1 <= position <= len(references[chrom]):
                    raise ValueError(f"Mpileup position is outside reference: {chrom}:{position}")
                if position != observed[chrom] + 1:
                    raise ValueError(f"Nonsequential/duplicate mpileup position: {chrom}:{position}")
                expected_ref = references[chrom][position - 1]
                if ref.upper() != expected_ref:
                    raise ValueError(
                        f"Mpileup/reference disagreement at {chrom}:{position}: {ref}!={expected_ref}"
                    )
                evidence = parse_read_bases(expected_ref, read_bases, depth)
                writer.writerow({
                    "Source": args.label,
                    "Reference": chrom,
                    "Position": position,
                    "RefBase": expected_ref,
                    "PileupDepth": depth,
                    **evidence,
                    "SNVNonRefFraction": f"{float(evidence['SNVNonRefFraction']):.8f}",
                    "MinMapQ": args.min_mapq,
                    "MinBaseQ": args.min_baseq,
                    "ReferenceSHA256": ref_sha,
                    "BAM": str(args.bam.resolve()),
                })
                observed[chrom] += 1
        return_code = process.wait()
        errors.seek(0)
        stderr = errors.read()
        if return_code != 0:
            raise RuntimeError(f"samtools mpileup failed ({return_code}): {stderr.strip()}")
        missing = {
            chrom: len(seq) - observed[chrom]
            for chrom, seq in references.items()
            if observed[chrom] != len(seq)
        }
        if missing:
            raise RuntimeError(f"mpileup -aa did not emit every reference position: {missing}")
        tmp.replace(args.output)
        print(
            f"[DONE] {args.output} source={args.label} mapq>={args.min_mapq} "
            f"baseq>={args.min_baseq}",
            flush=True,
        )
        return 0
    except Exception:
        process.kill()
        process.wait()
        tmp.unlink(missing_ok=True)
        raise
    finally:
        process.stdout.close()
        errors.close()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
