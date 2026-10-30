#!/usr/bin/env python3
"""Add base-resolved gene evidence without modifying existing analysis results.

Shared by Step3, Step4 and backfill. Sidecars record reference sequence SHA256,
source identity, generation, options, tool version and output SHA256.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from rrmg_pileup_base_counts import reference_inventory

SCHEMA = 1


def zero_target_design(row):
    total = row.get("total_pairs") or row.get("computed_total_pairs")
    abundance = row.get("target_abundance")
    return total not in (None, "") and abundance not in (None, "") and round(float(total) * float(abundance)) == 0


def validated_evidence(path, reference_sha, reference_name):
    """Validate a derived file without needing the large source FASTA/BAM."""
    path = Path(path)
    saved = json.loads(marker_path(path).read_text(encoding="utf-8"))
    signature = saved["signature"]
    if signature["schema"] != SCHEMA or signature["reference_sha256"] != reference_sha or signature["reference_name"] != reference_name:
        raise ValueError(f"Evidence/reference fingerprint mismatch: {path}")
    if saved["output_sha256"] != sha256(path):
        raise ValueError(f"Evidence checksum mismatch: {path}")
    source = signature.get("source")
    if source and Path(source["path"]).is_file() and identity(source["path"]) != source:
        raise ValueError(f"Evidence source changed; backfill again: {path}")
    return signature


def validate_base_table(path, name, sequence, signature):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        n = 0
        for n, row in enumerate(rows, 1):
            if n > len(sequence) or row["Reference"] != name or int(row["Position"]) != n or row["RefBase"] != sequence[n - 1]:
                raise ValueError(f"Base-count reference/coordinate mismatch: {path}")
            if row["ReferenceSHA256"] != signature["reference_sha256"] or int(row["MinMapQ"]) != signature["min_mapq"] or int(row["MinBaseQ"]) != signature["min_baseq"]:
                raise ValueError(f"Base-count parameter mismatch: {path}")
            if sum(int(row[b]) for b in "ACGTN") != int(row["BaseCountDepth"]):
                raise ValueError(f"Base-count depth mismatch: {path}")
        if n != len(sequence):
            raise ValueError(f"Incomplete per-base evidence: {path}")


def interval_length(intervals):
    end = total = 0
    for a, b in sorted(intervals):
        total += max(0, b - max(a, end))
        end = max(end, b)
    return total


def contig_details(path, sequence):
    """Recovery classification and exact changes; PAF coordinates are 0-based."""
    by_contig, changes = {}, []
    tokens = re.compile(r"=([A-Za-z]+)|:(\d+)|\*([A-Za-z])([A-Za-z])|\+([A-Za-z]+)|-([A-Za-z]+)")
    with Path(path).open() as stream:
        for raw in stream:
            f = raw.rstrip().split("\t")
            if len(f) < 12:
                raise ValueError(f"Malformed cs PAF: {path}")
            start, end, block, mq = int(f[7]), int(f[8]), int(f[10]), int(f[11])
            ident = 100 * int(f[9]) / block if block else 0
            passed = mq >= 20 and block >= 150 and ident >= 95
            if passed:
                by_contig.setdefault(f[0], []).append((start, end))
            cs = next((v[5:] for v in f[12:] if v.startswith("cs:Z:")), None)
            if cs is None:
                raise ValueError(f"cs tag is absent: {path}")
            cursor, consumed, qconsumed = start, 0, 0
            for match in tokens.finditer(cs):
                if match.start() != consumed:
                    raise ValueError(f"Unsupported cs operation: {cs}")
                consumed = match.end()
                exact, length, ref, alt, ins, deletion = match.groups()
                if exact or length:
                    n = len(exact) if exact else int(length)
                    if exact and sequence[cursor:cursor + n] != exact.upper():
                        raise ValueError("cs exact-match sequence disagrees with reference")
                    cursor += n; qconsumed += n
                    continue
                event = "SNV" if ref else "INS" if ins else "DEL"
                ref_allele = ref.upper() if ref else deletion.upper() if deletion else "-"
                alt_allele = alt.upper() if alt else ins.upper() if ins else "-"
                position = cursor + 1 if event != "INS" else cursor
                if event != "INS" and sequence[cursor:cursor + len(ref_allele)] != ref_allele:
                    raise ValueError("cs difference disagrees with reference")
                changes.append(dict(ContigID=f[0], GenePosition=position,
                                    EventType=event, RefAllele=ref_allele, QueryAllele=alt_allele,
                                    Strand=f[4], MapQ=mq, Identity=f"{ident:.6f}",
                                    PassRecoveryFilter="Y" if passed else "N"))
                cursor += len(ref_allele) if event != "INS" else 0
                qconsumed += len(alt_allele) if event != "DEL" else 0
            if consumed != len(cs) or cursor != end or qconsumed != int(f[3]) - int(f[2]):
                raise ValueError("cs length does not match PAF intervals")
    breadths = [100 * interval_length(v) / len(sequence) for v in by_contig.values()]
    maximum = max(breadths, default=0)
    union = 100 * interval_length([iv for v in by_contig.values() for iv in v]) / len(sequence)
    n = len(by_contig)
    category = ("NOT_DETECTED" if not n else "FULL_SINGLE" if maximum >= 90 and n == 1
                else "MULTI_CONTIG_WITH_FULL" if maximum >= 90
                else "FRAGMENTED_COMPLETE" if union >= 90 and n >= 2 else "FRAGMENTED_PARTIAL")
    return dict(RecoveryClass=category, QualifyingContigN=n, MaxSingleContigBreadth=f"{maximum:.6f}",
                CombinedContigBreadth=f"{union:.6f}", DifferenceEventN=len(changes)), changes


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def version(tool):
    r = subprocess.run([tool, "--version"], capture_output=True, text=True, check=True, errors="replace")
    return r.stdout.splitlines()[0] if r.stdout.splitlines() else "unreported"


def marker_path(output):
    return Path(str(output) + ".evidence.json")


def is_current(output, signature):
    try:
        saved = json.loads(marker_path(output).read_text(encoding="utf-8"))
        return saved["signature"] == signature and saved["output_sha256"] == sha256(output)
    except (OSError, ValueError, KeyError):
        return False


def mark(output, signature):
    path = marker_path(output)
    tmp = Path(str(path) + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(dict(signature=signature, output_sha256=sha256(output)), indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def ensure_contigs(reference, contigs, output, generation="", threads=4, minimap2="minimap2", force=False):
    reference, contigs, output = map(Path, (reference, contigs, output))
    if not contigs.is_file() or not contigs.stat().st_size:
        raise ValueError(f"Missing/empty contigs: {contigs}")
    refs, ref_sha = reference_inventory(reference)
    signature = dict(schema=SCHEMA, kind="contig_cs", reference_sha256=ref_sha,
                     reference_name=next(iter(refs)), source=identity(contigs), generation=generation,
                     preset="asm5", cs="long", cigar=True, tool=version(minimap2))
    if not force and is_current(output, signature):
        return "SKIP_CURRENT"
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(output) + f".tmp.{os.getpid()}")
    try:
        with tmp.open("w", encoding="utf-8", newline="") as stream, Path(str(output) + ".minimap2.log").open("w") as log:
            subprocess.run([minimap2, "-t", str(threads), "-x", "asm5", "-c", "--cs=long", str(reference), str(contigs)], stdout=stream, stderr=log, check=True)
        # Empty PAF is an explicit no-hit, with the sidecar proving command success.
        with tmp.open() as stream:
            for line in stream:
                fields = line.rstrip().split("\t")
                if len(fields) < 12 or not any(f.startswith("cs:Z:") for f in fields[12:]) or not any(f.startswith("cg:Z:") for f in fields[12:]):
                    raise ValueError("minimap2 output lacks base-level cs/cg tags")
        tmp.replace(output)
        mark(output, signature)
        return "OK"
    finally:
        tmp.unlink(missing_ok=True)


def ensure_reads(reference, prefix, label, generation="", samtools="samtools", force=False, empty=False):
    reference, prefix = Path(reference), Path(prefix)
    bam = Path(str(prefix) + ".bam")
    if not empty and (not bam.is_file() or not bam.stat().st_size):
        raise ValueError(f"Missing/empty BAM; cannot infer zero coverage: {bam}")
    refs, ref_sha = reference_inventory(reference)
    tool_version = version(samtools) if not empty else "not_run_zero_design"
    source = None if empty else identity(bam)
    states = []
    for suffix, mq, bq in (("raw", 0, 0), ("q20", 20, 20)):
        output = Path(str(prefix) + f".base_counts.{suffix}.tsv.gz")
        signature = dict(schema=SCHEMA, kind="read_base_counts", reference_sha256=ref_sha,
                         reference_name=next(iter(refs)), source=source, generation=generation,
                         min_mapq=mq, min_baseq=bq, baq=False, count_overlap=True,
                         count_orphans=True, exclude_flags="0xF04", depth_cap=0,
                         label=label, tool=tool_version, explicit_zero=empty,
                         parser_sha256=sha256(Path(__file__).with_name("rrmg_pileup_base_counts.py")))
        if not force and is_current(output, signature):
            states.append("SKIP_CURRENT")
            continue
        cmd = [sys.executable, str(Path(__file__).with_name("rrmg_pileup_base_counts.py")),
               "--reference", str(reference), "--output", str(output), "--label", label,
               "--samtools", samtools, "--min-mapq", str(mq), "--min-baseq", str(bq)]
        cmd += ["--empty"] if empty else ["--bam", str(bam)]
        subprocess.run(cmd, check=True, stdout=sys.stderr)
        mark(output, signature)
        states.append("OK")
    return "/".join(states)


def generation_status(current, marker):
    if not current or not marker:
        return "UNVERIFIED"
    return "MATCH" if current == marker else "MISMATCH"


def text_marker(path):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return ""


def mix_generation(run):
    try:
        with (run / "design/read_mix_audit.tsv").open() as stream:
            return next((r[1] for r in csv.reader(stream, delimiter="\t") if len(r) >= 2 and r[0] == "mix_generation_id"), "")
    except OSError:
        return ""


def select_plan(plan, spec):
    with Path(plan).open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise ValueError("Empty plan")
    ids, names = set(), set()
    for row in rows:
        index, name = int(row["run_index"]), row["run_id"]
        if index <= 0 or index in ids or name in names or not name or "/" in name or "\\" in name or name in (".", ".."):
            raise ValueError("Duplicate/invalid plan ID or unsafe run_id")
        ids.add(index); names.add(name)
    if spec == "all":
        return rows
    selected = set()
    for part in spec.split(","):
        if not re.fullmatch(r"[1-9]\d*(?:-[1-9]\d*)?", part):
            raise ValueError("Invalid --rows")
        ends = list(map(int, part.split("-")))
        if ends[0] > ends[-1]:
            raise ValueError("Descending --rows")
        selected.update(range(ends[0], ends[-1] + 1))
    if selected - ids:
        raise ValueError(f"Indices absent from plan: {sorted(selected - ids)}")
    return [row for row in rows if int(row["run_index"]) in selected]


def backfill(args):
    wd = args.wd.resolve()
    reference = args.reference or wd / "ref/CYP51A_wt.fasta"
    rows = select_plan(args.plan or wd / "input/rrmg/rrmg_permutation_plan.tsv", args.rows)
    # Do not write ad-hoc diagnostics inside the collector-managed figure_data.
    status = args.status or wd / "result/CYP51A_evidence_backfill_status.tsv"
    status.parent.mkdir(parents=True, exist_ok=True)
    fields = ["run_index", "run_id", "contig_cs", "metagenomic_counts", "synthetic_counts",
              "contig_generation", "metagenomic_generation", "synthetic_generation", "error"]
    failures = 0
    with status.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            run = wd / "temp/rrmg_permutation" / row["run_id"]
            current = mix_generation(run)
            asm_gen = text_marker(run / "asm/RRMG_MIX_GENERATION")
            mix_gen = text_marker(run / "readmap/reads_to_GENE.mix_generation.txt")
            target_gen = text_marker(run / "readmap/target_only_reads_to_GENE.mix_generation.txt")
            record = dict(run_index=row["run_index"], run_id=row["run_id"],
                          contig_generation=generation_status(current, asm_gen),
                          metagenomic_generation=generation_status(current, mix_gen),
                          synthetic_generation=generation_status(current, target_gen))
            errors = []
            print(f"[ROW] {row['run_index']} {row['run_id']}", flush=True)
            jobs = [("contig_cs", lambda: ensure_contigs(reference, run / "asm/contigs.fasta", run / "asm/contigs.fasta.to_GENE.cs.paf", asm_gen, args.threads, args.minimap2, args.force)),
                    ("metagenomic_counts", lambda: ensure_reads(reference, run / "readmap/reads_to_GENE", "metagenomic", mix_gen, args.samtools, args.force)),
                    ("synthetic_counts", lambda: ensure_reads(reference, run / "readmap/target_only_reads_to_GENE", "synthetic_target", target_gen, args.samtools, args.force,
                        empty=zero_target_design(row)))]
            for name, job in jobs:
                try:
                    record[name] = job()
                except (OSError, ValueError, subprocess.SubprocessError) as exc:
                    record[name] = "FAILED"
                    errors.append(f"{name}: {exc}")
            record["error"] = "; ".join(errors)
            failures += bool(errors)
            writer.writerow(record); stream.flush()
    print(f"[DONE] rows={len(rows)} failed_rows={failures} status={status}", flush=True)
    return int(failures > 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("contigs")
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument("--contigs", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--generation", default="")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--force", action="store_true")
    p = sub.add_parser("reads")
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument("--prefix", required=True, type=Path)
    p.add_argument("--label", required=True)
    p.add_argument("--generation", default="")
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--empty", action="store_true")
    p.add_argument("--force", action="store_true")
    p = sub.add_parser("backfill")
    p.add_argument("--wd", "-w", required=True, type=Path)
    p.add_argument("--code", "-c", type=Path, help="compatibility option; helpers are loaded beside this script")
    p.add_argument("--reference", type=Path)
    p.add_argument("--plan", type=Path)
    p.add_argument("--rows", "-r", default="all")
    p.add_argument("--threads", "-t", type=int, default=4)
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--status", type=Path)
    p.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if hasattr(args, "threads") and args.threads < 1:
        parser.error("threads must be positive")
    if args.command == "contigs":
        print(ensure_contigs(args.reference, args.contigs, args.output, args.generation, args.threads, args.minimap2, args.force))
    elif args.command == "reads":
        print(ensure_reads(args.reference, args.prefix, args.label, args.generation, args.samtools, args.force, args.empty))
    else:
        return backfill(args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
