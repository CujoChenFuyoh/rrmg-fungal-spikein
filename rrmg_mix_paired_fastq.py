#!/usr/bin/env python3
"""Count, exact-subsample, and reproducibly interleave paired FASTQ files.

The ``mix`` command always writes one gzip member per mate.  It therefore also
normalizes BGZF or concatenated-gzip inputs into a stream that readers such as
fastp can consume without stopping at an internal EOF marker.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import uuid


class FastqError(ValueError):
    """Invalid, truncated, or unexpectedly sized FASTQ input."""


def open_fastq(path: Path):
    return gzip.open(path, "rb") if path.name.endswith(".gz") else path.open("rb")


def fastq_records(handle, path: Path):
    record_number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        record_number += 1
        if not sequence or not plus or not quality:
            raise FastqError(f"Truncated FASTQ record {record_number}: {path}")
        if not header.startswith(b"@") or not plus.startswith(b"+"):
            raise FastqError(f"Invalid FASTQ record {record_number}: {path}")
        if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
            raise FastqError(f"Sequence/quality length mismatch in record {record_number}: {path}")
        yield header, sequence, plus, quality


def read_key(record) -> bytes:
    """Return the mate-independent first token of a FASTQ header."""
    header = record[0].rstrip(b"\r\n").split(None, 1)[0]
    if header.startswith(b"@"):
        header = header[1:]
    if header.endswith((b"/1", b"/2")):
        header = header[:-2]
    return header


def paired_records(r1: Path, r2: Path):
    with open_fastq(r1) as handle1, open_fastq(r2) as handle2:
        records1 = fastq_records(handle1, r1)
        records2 = fastq_records(handle2, r2)
        pair_number = 0
        while True:
            record1 = next(records1, None)
            record2 = next(records2, None)
            if record1 is None and record2 is None:
                return
            pair_number += 1
            if record1 is None or record2 is None:
                raise FastqError(
                    f"Mate files contain different numbers of records near pair {pair_number}: {r1} / {r2}"
                )
            if read_key(record1) != read_key(record2):
                raise FastqError(
                    f"Mate names differ at pair {pair_number}: "
                    f"{record1[0].decode(errors='replace').strip()} / "
                    f"{record2[0].decode(errors='replace').strip()}"
                )
            yield record1, record2


def count_pairs(r1: Path, r2: Path) -> int:
    return sum(1 for _ in paired_records(r1, r2))


def selected_pairs(records, source_total: int, wanted: int, rng: random.Random, label: str):
    """Uniformly choose ``wanted`` records from a stream of known exact size."""
    if source_total < 0 or wanted < 0 or wanted > source_total:
        raise FastqError(f"Invalid {label} counts: source={source_total}, requested={wanted}")
    remaining_source = source_total
    remaining_wanted = wanted
    for pair_index in range(1, source_total + 1):
        record = next(records, None)
        if record is None:
            observed = pair_index - 1
            raise FastqError(
                f"{label} FASTQ ended early: declared {source_total} pairs, observed {observed}"
            )
        if remaining_wanted and rng.randrange(remaining_source) < remaining_wanted:
            remaining_wanted -= 1
            yield record
        remaining_source -= 1
    if remaining_wanted:
        raise FastqError(f"Could not select {wanted} {label} pairs")
    if next(records, None) is not None:
        raise FastqError(f"{label} FASTQ has more than the declared {source_total} pairs")


class CompressedWriter:
    def __init__(self, destination: Path, compressor: str, threads: int, level: int):
        self.destination = destination
        self.compressor = compressor
        self.threads = threads
        self.level = level
        self.temp_path: Path | None = None
        self.raw = None
        self.stream = None
        self.process = None
        self.stderr = b""

    def open(self):
        self.destination.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(prefix="." + self.destination.name + ".", suffix=".tmp", dir=self.destination.parent)
        os.close(fd)
        self.temp_path = Path(name)
        self.raw = self.temp_path.open("wb")
        if self.compressor == "pigz":
            self.process = subprocess.Popen(
                ["pigz", f"-{self.level}", "-n", "-p", str(self.threads), "-c"],
                stdin=subprocess.PIPE,
                stdout=self.raw,
                stderr=subprocess.PIPE,
            )
            self.stream = self.process.stdin
        else:
            self.stream = gzip.GzipFile(
                filename="", mode="wb", compresslevel=self.level, fileobj=self.raw, mtime=0
            )
        return self

    def write(self, record):
        for line in record:
            self.stream.write(line)

    def close(self):
        if self.stream is not None:
            self.stream.close()
            self.stream = None
        if self.process is not None:
            self.stderr = self.process.stderr.read()
            returncode = self.process.wait()
            self.process.stderr.close()
            self.process = None
            if returncode:
                raise FastqError("pigz failed: " + self.stderr.decode(errors="replace").strip())
        if self.raw is not None:
            self.raw.close()
            self.raw = None

    def abort(self):
        try:
            if self.stream is not None:
                self.stream.close()
        except (BrokenPipeError, OSError):
            pass
        if self.process is not None:
            self.process.kill()
            self.process.wait()
        if self.raw is not None:
            self.raw.close()
        if self.temp_path is not None:
            self.temp_path.unlink(missing_ok=True)


def choose_compressor(requested: str) -> str:
    if requested == "auto":
        return "pigz" if shutil.which("pigz") else "python"
    if requested == "pigz" and not shutil.which("pigz"):
        raise FastqError("--compressor pigz was requested, but pigz is not on PATH")
    return requested


def ensure_distinct(paths):
    resolved = [path.resolve() for path in paths if path is not None]
    if len(resolved) != len(set(resolved)):
        raise FastqError("Input and output FASTQ paths must all be distinct")


def write_audit(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix="." + path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["key", "value"])
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def integer(value, label: str) -> int:
    try:
        number = int(round(float(value)))
    except (TypeError, ValueError) as error:
        raise FastqError(f"Invalid {label}: {value}") from error
    if number < 0:
        raise FastqError(f"Invalid {label}: {value}")
    return number


def validate_mix_audit(
    path: Path,
    expected_pairs: int,
    r1: Path | None = None,
    r2: Path | None = None,
    expected_background_pairs: int | None = None,
    expected_target_pairs: int | None = None,
):
    if not path.is_file():
        raise FastqError(f"Read-mix audit not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    values = {row.get("key", ""): row.get("value", "") for row in rows}
    observed = integer(values.get("total_pairs_written"), "audited total_pairs_written")
    background = integer(values.get("background_pairs_written"), "audited background_pairs_written")
    target = integer(values.get("target_pairs_written"), "audited target_pairs_written")
    members = integer(values.get("gzip_members_per_output"), "audited gzip_members_per_output")
    generation = values.get("mix_generation_id", "")
    if len(generation) != 32 or any(character not in "0123456789abcdef" for character in generation):
        raise FastqError(f"Invalid or missing mix_generation_id in read-mix audit: {path}")
    if observed != expected_pairs or background + target != observed or members != 1:
        raise FastqError(
            f"Read-mix audit mismatch: expected_pairs={expected_pairs}, "
            f"audited_pairs={observed}, background_plus_target={background + target}, "
            f"gzip_members={members}"
        )
    if expected_background_pairs is not None and background != expected_background_pairs:
        raise FastqError(
            f"Audited background pairs differ from design: {background} / {expected_background_pairs}"
        )
    if expected_target_pairs is not None and target != expected_target_pairs:
        raise FastqError(f"Audited target pairs differ from design: {target} / {expected_target_pairs}")
    if (r1 is None) != (r2 is None):
        raise FastqError("Both --r1 and --r2 are required when validating mixed FASTQ files")
    if r1 is not None and r2 is not None:
        for label, supplied, key_path, key_bytes in (
            ("R1", r1.resolve(), "out_r1", "out_r1_bytes"),
            ("R2", r2.resolve(), "out_r2", "out_r2_bytes"),
        ):
            audited = Path(values.get(key_path, "")).resolve()
            if supplied != audited:
                raise FastqError(f"Mixed {label} path differs from its audit: {supplied} / {audited}")
            if not supplied.is_file():
                raise FastqError(f"Mixed {label} is missing: {supplied}")
            audited_bytes = integer(values.get(key_bytes), f"audited {key_bytes}")
            if supplied.stat().st_size != audited_bytes:
                raise FastqError(
                    f"Mixed {label} size differs from its audit: "
                    f"{supplied.stat().st_size} / {audited_bytes} bytes"
                )
    print(
        f"[MIX-AUDIT] expected_pairs={expected_pairs} audited_pairs={observed} "
        f"gzip_members={members} generation={generation}"
    )


def validate_fastp(path: Path, expected_pairs: int):
    if not path.is_file():
        raise FastqError(f"fastp JSON not found: {path}")
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
        observed_reads = integer(report["summary"]["before_filtering"]["total_reads"], "fastp total_reads")
    except (KeyError, json.JSONDecodeError) as error:
        raise FastqError(f"Invalid fastp JSON: {path}: {error}") from error
    expected_reads = 2 * expected_pairs
    if observed_reads != expected_reads:
        raise FastqError(
            f"fastp did not consume the complete paired input: expected {expected_reads} reads "
            f"({expected_pairs} pairs), observed {observed_reads}. Do not use this QC/assembly."
        )
    print(f"[FASTP-AUDIT] expected_reads={expected_reads} observed_reads={observed_reads}")


def mix(args):
    background_r1 = args.background_r1.resolve()
    background_r2 = args.background_r2.resolve()
    target_r1 = args.target_r1.resolve() if args.target_r1 else None
    target_r2 = args.target_r2.resolve() if args.target_r2 else None
    output_r1 = args.out_r1.resolve()
    output_r2 = args.out_r2.resolve()
    audit = args.audit.resolve()

    if args.target_pairs and (target_r1 is None or target_r2 is None):
        raise FastqError("Target FASTQs are required when --target-pairs is positive")
    if not args.target_pairs and (target_r1 is not None or target_r2 is not None):
        raise FastqError("Do not provide target FASTQs when --target-pairs is zero")
    if args.target_source_pairs is None:
        args.target_source_pairs = args.target_pairs
    ensure_distinct([background_r1, background_r2, target_r1, target_r2, output_r1, output_r2])

    compressor = choose_compressor(args.compressor)
    background = selected_pairs(
        iter(paired_records(background_r1, background_r2)),
        args.background_source_pairs,
        args.background_pairs,
        random.Random(args.seed + 101),
        "background",
    )
    target = selected_pairs(
        iter(paired_records(target_r1, target_r2)) if args.target_pairs else iter(()),
        args.target_source_pairs,
        args.target_pairs,
        random.Random(args.seed + 202),
        "target",
    )

    writer1 = CompressedWriter(output_r1, compressor, args.threads, args.compression_level).open()
    writer2 = CompressedWriter(output_r2, compressor, args.threads, args.compression_level).open()
    background_remaining = args.background_pairs
    target_remaining = args.target_pairs
    source_switches = 0
    previous_source = None
    first_target = None
    last_target = None
    interleave_rng = random.Random(args.seed + 303)
    total_pairs = background_remaining + target_remaining
    try:
        for output_index in range(1, total_pairs + 1):
            remaining = background_remaining + target_remaining
            use_target = bool(target_remaining and interleave_rng.randrange(remaining) < target_remaining)
            if use_target:
                pair = next(target)
                target_remaining -= 1
                source = "target"
                first_target = output_index if first_target is None else first_target
                last_target = output_index
            else:
                pair = next(background)
                background_remaining -= 1
                source = "background"
            if previous_source is not None and source != previous_source:
                source_switches += 1
            previous_source = source
            writer1.write(pair[0])
            writer2.write(pair[1])

        # Resume each selector after its final yield so it validates the full
        # declared source size and paired EOF state.
        if next(background, None) is not None or next(target, None) is not None:
            raise FastqError("Internal selection error: unexpected extra selected pair")
        writer1.close()
        writer2.close()
        os.replace(writer1.temp_path, output_r1)
        os.replace(writer2.temp_path, output_r2)
    except Exception:
        writer1.abort()
        writer2.abort()
        raise

    rows = [
        ("mix_generation_id", uuid.uuid4().hex),
        ("method", "uniform_exact_subsample_then_uniform_paired_interleave"),
        ("background_r1", background_r1),
        ("background_r2", background_r2),
        ("background_source_pairs", args.background_source_pairs),
        ("background_pairs_written", args.background_pairs),
        ("target_r1", target_r1 or ""),
        ("target_r2", target_r2 or ""),
        ("target_source_pairs", args.target_source_pairs),
        ("target_pairs_written", args.target_pairs),
        ("total_pairs_written", total_pairs),
        ("mix_seed", args.seed),
        ("source_switches", source_switches),
        ("first_target_pair_index", first_target if first_target is not None else "NA"),
        ("last_target_pair_index", last_target if last_target is not None else "NA"),
        ("compressor", compressor),
        ("gzip_members_per_output", 1),
        ("out_r1", output_r1),
        ("out_r2", output_r2),
        ("out_r1_bytes", output_r1.stat().st_size),
        ("out_r2_bytes", output_r2.stat().st_size),
    ]
    write_audit(audit, rows)
    print(
        f"[MIX] total_pairs={total_pairs} background_pairs={args.background_pairs} "
        f"target_pairs={args.target_pairs} switches={source_switches} compressor={compressor}"
    )


def parser():
    top = argparse.ArgumentParser(description=__doc__)
    commands = top.add_subparsers(dest="command", required=True)
    count = commands.add_parser("count", help="Validate paired FASTQ and print its pair count")
    count.add_argument("--r1", type=Path, required=True)
    count.add_argument("--r2", type=Path, required=True)

    command = commands.add_parser("mix", help="Write an exact, paired, randomly interleaved single-member gzip")
    command.add_argument("--background-r1", type=Path, required=True)
    command.add_argument("--background-r2", type=Path, required=True)
    command.add_argument("--background-source-pairs", type=int, required=True)
    command.add_argument("--background-pairs", type=int, required=True)
    command.add_argument("--target-r1", type=Path)
    command.add_argument("--target-r2", type=Path)
    command.add_argument("--target-source-pairs", type=int)
    command.add_argument("--target-pairs", type=int, default=0)
    command.add_argument("--out-r1", type=Path, required=True)
    command.add_argument("--out-r2", type=Path, required=True)
    command.add_argument("--audit", type=Path, required=True)
    command.add_argument("--seed", type=int, required=True)
    command.add_argument("--threads", type=int, default=1)
    command.add_argument("--compression-level", type=int, choices=range(1, 10), default=4)
    command.add_argument("--compressor", choices=("auto", "pigz", "python"), default="auto")

    command = commands.add_parser("validate-mix", help="Require an exact, single-member Step 2 audit")
    command.add_argument("--audit", type=Path, required=True)
    command.add_argument("--expected-pairs", required=True)
    command.add_argument("--r1", type=Path)
    command.add_argument("--r2", type=Path)
    command.add_argument("--expected-background-pairs")
    command.add_argument("--expected-target-pairs")

    command = commands.add_parser("validate-fastp", help="Require fastp to have consumed every planned input pair")
    command.add_argument("--json", dest="json_path", type=Path, required=True)
    command.add_argument("--expected-pairs", required=True)
    return top


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == "count":
            print(count_pairs(args.r1.resolve(), args.r2.resolve()))
        elif args.command == "mix":
            mix(args)
        elif args.command == "validate-mix":
            validate_mix_audit(
                args.audit.resolve(),
                integer(args.expected_pairs, "expected pairs"),
                args.r1,
                args.r2,
                integer(args.expected_background_pairs, "expected background pairs")
                if args.expected_background_pairs is not None
                else None,
                integer(args.expected_target_pairs, "expected target pairs")
                if args.expected_target_pairs is not None
                else None,
            )
        else:
            validate_fastp(args.json_path.resolve(), integer(args.expected_pairs, "expected pairs"))
        return 0
    except (FastqError, OSError, subprocess.SubprocessError) as error:
        parser().exit(1, "[ERROR] " + str(error) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
