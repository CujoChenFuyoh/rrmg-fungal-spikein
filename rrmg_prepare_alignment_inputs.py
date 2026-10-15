#!/usr/bin/env python3
"""Prepare portable gene-alignment figures, without changing Step1-5 results.

Remaps existing contigs, QC reads and retained target-only reads to ONE gene
reference. No read regeneration, reference editing or mutation calling occurs.
Site annotations are validated against the supplied catalogue and reference.
"""
import argparse
from contextlib import contextmanager
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path


def input_file(value):
    """Reject unset shell variables before Path('') can become the current directory."""
    if not value.strip():
        raise argparse.ArgumentTypeError("Empty file path; set SITES_SOURCE/CATALOGUE_SOURCE (runbook section 0 or 14).")
    path = Path(value)
    if not path.is_file() or not path.stat().st_size:
        raise argparse.ArgumentTypeError(f"Expected a non-empty file, not a directory or missing path: {value}")
    return path


OUTPUT_MARKER = ".rrmg_alignment_outputs.json"


def output_inventory(path):
    result = {}
    for item in sorted(path.rglob("*")):
        if item.is_symlink():
            raise ValueError("Output symlinks are not supported: " + str(item))
        if item.is_file() and item.name != OUTPUT_MARKER:
            digest = hashlib.sha256()
            with item.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            result[item.relative_to(path).as_posix()] = digest.hexdigest()
    return result


def check_output(path):
    if path.is_symlink():
        raise ValueError("Output directory cannot be a symlink: " + str(path))
    if not path.exists():
        return
    marker = path / OUTPUT_MARKER
    if not marker.is_file():
        raise ValueError("Existing output is unmanaged; preserve/move it aside first: " + str(path))
    saved = json.loads(marker.read_text(encoding="utf-8"))
    if saved.get("kind") != "rrmg_alignment_inputs" or saved.get("sha256") != output_inventory(path):
        raise ValueError("Output contains added, changed or missing files; preserve/move it aside first: " + str(path))


@contextmanager
def managed_output(path):
    """Build beside a fixed destination, then replace only verified tool-owned outputs."""
    path = Path(path).absolute()
    path.parent.mkdir(parents=True, exist_ok=True)
    lock = path.with_name(path.name + ".lock")
    try:
        lock.mkdir()
    except FileExistsError:
        raise ValueError(f"Output is locked: {lock}; check for a running process before removing the empty lock directory.") from None
    try:
        check_output(path)
        with tempfile.TemporaryDirectory(prefix="." + path.name + "-stage-", dir=path.parent) as tmp:
            stage = Path(tmp) / "output"
            stage.mkdir()
            yield stage
            (stage / OUTPUT_MARKER).write_text(json.dumps(dict(kind="rrmg_alignment_inputs", sha256=output_inventory(stage)), indent=2) + "\n", encoding="utf-8")
            check_output(path)
            # Keep the previous results OUTSIDE auto-cleaned staging so a failed
            # rollback can never delete the last successful run.
            backup = path.with_name("." + path.name + "-previous-" + uuid.uuid4().hex)
            had_old = path.exists()
            if had_old:
                path.rename(backup)
            try:
                stage.rename(path)
            except OSError:
                if had_old:
                    try:
                        backup.rename(path)
                    except OSError as exc:
                        raise RuntimeError("Previous results preserved at: " + str(backup)) from exc
                raise
            if had_old:
                if backup.resolve().parent != path.parent.resolve():
                    raise RuntimeError("Unsafe backup cleanup path: " + str(backup))
                shutil.rmtree(backup)
    finally:
        lock.rmdir()


def table(path):
    if not Path(path).is_file():
        raise ValueError(f"Expected a TSV file, not a directory or missing path: {path}")
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows, fields):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    headers = [line[1:] for line in lines if line.startswith(">")]
    if len(headers) != 1:
        raise ValueError("Alignment figures require exactly one reference sequence.")
    sequence = "".join(line.strip() for line in lines if not line.startswith(">" )).upper()
    if not sequence or re.search("[^ACGTRYSWKMBDHVN]", sequence):
        raise ValueError("Empty/invalid reference FASTA.")
    return headers[0], sequence


def validated_sites(gene, site_file, catalogue):
    """Use exact chromosome codon coordinates, not amino-acid-position * 3.

    The NCBI interval in the FASTA header establishes orientation. Catalogue
    status and reference codon must agree. Malformed/inconsistent sites are
    excluded with an explicit reason; coordinates are never reconstructed.
    """
    header, seq = fasta(gene)
    match = re.search(r"(NC_\d+\.\d+):(c?)(\d+)-(\d+)", header)
    if not match:
        raise ValueError("Site validation needs an NCBI chromosome interval in the target FASTA header.")
    chrom, reverse, a, b = match.groups()
    a, b = int(a), int(b)
    start, end = sorted((a, b))
    strand = "-" if reverse else "+"
    if len(seq) != end - start + 1 or (strand == "-" and a < b) or (strand == "+" and a > b):
        raise ValueError("Reference interval/orientation/sequence length disagree.")
    catalogue_rows = table(catalogue)
    by_label = {}
    for row in catalogue_rows:
        if row.get("gene", "").upper() == "CYP51A":
            by_label.setdefault(row.get("aa_change", ""), []).append(row)
    clean, audit = [], []
    seen = set()
    for site in table(site_file):
        label = site.get("label", "")
        errors, accepted = [], None
        for row in by_label.get(label, []):
            try:
                if row.get("status", "").lower() != "ok":
                    raise ValueError("catalogue_status=" + row.get("status", "missing"))
                if row.get("chrom") != chrom or row.get("strand") != strand:
                    raise ValueError("chromosome/strand mismatch")
                raw = row.get("codon_positions", "")
                if not re.fullmatch(r"\d+\|\d+\|\d+", raw):
                    raise ValueError("codon_positions must contain three exact pipe-separated coordinates")
                coords = [int(x) for x in raw.split("|")]
                positions = [end - p + 1 if reverse else p - start + 1 for p in coords]
                if len(set(positions)) != 3 or not all(1 <= p <= len(seq) for p in positions):
                    raise ValueError("codon outside gene interval")
                observed = "".join(seq[p - 1] for p in positions)
                if observed != row.get("ref_codon", "").upper():
                    raise ValueError("reference codon mismatch")
                if (float(site["plot_start"]), float(site["plot_end"])) != (min(positions), max(positions)):
                    raise ValueError("existing plot coordinates disagree with catalogue")
                if "plot_mid" in site and abs(float(site["plot_mid"]) - (min(positions) + max(positions)) / 2) > 0.001:
                    raise ValueError("existing midpoint disagrees with catalogue")
                # The legacy interval-based plot cannot represent a codon split by an intron.
                if max(positions) - min(positions) != 2:
                    raise ValueError("split codon needs a non-contiguous annotation renderer")
                accepted = dict(label=label, gene="CYP51A", chrom=chrom, chr_start=min(coords), chr_end=max(coords),
                                plot_start=min(positions), plot_end=max(positions),
                                plot_mid=(min(positions) + max(positions)) / 2,
                                mutation_type=site.get("mutation_type", "AA_substitution"),
                                is_key_site=site.get("is_key_site", "N"),
                                marker_id=row["marker_id"], ref_codon=observed,
                                reference_sha256=hashlib.sha256(seq.encode()).hexdigest())
                break
            except (ValueError, KeyError) as exc:
                errors.append(str(exc))
        if accepted and label not in seen:
            clean.append(accepted)
            seen.add(label)
            audit.append(dict(label=label, status="included", reason="catalogue coordinates and reference codon verified"))
        else:
            audit.append(dict(label=label, status="excluded", reason="; ".join(sorted(set(errors))) or "no matching catalogue entry or duplicate"))
    if not clean:
        raise ValueError("No site annotations passed reference/catalogue validation.")
    return clean, audit


def paired_reads(folder, pattern):
    matches = sorted(folder.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"Expected one R1 matching {folder / pattern}; found {len(matches)}. Retain Step2 target reads and Step3 QC reads.")
    r1 = matches[0]
    r2 = r1.with_name(r1.name.replace("_1.fastq.gz", "_2.fastq.gz"))
    if not r2.is_file() or min(r1.stat().st_size, r2.stat().st_size) == 0:
        raise ValueError(f"Missing/empty paired FASTQ: {r1}, {r2}")
    return r1, r2


def command(cmd, out, log):
    with Path(out).open("wb") as handle, Path(log).open("wb") as err:
        subprocess.run([str(c) for c in cmd], stdout=handle, stderr=err, check=True)


def paf_stats(path, name, length):
    covered = set()
    count = 0
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12 or fields[5] != name or int(fields[6]) != length:
                raise ValueError(f"Invalid PAF/reference mismatch: {path}")
            start, end = int(fields[7]), int(fields[8])
            if not 0 <= start < end <= length:
                raise ValueError(f"PAF coordinate out of range: {path}")
            covered.update(range(start, end))
            count += 1
    return count, 100 * len(covered) / length


def map_reads(gene, reads, prefix, threads, name, length):
    paf = prefix.with_suffix(".paf")
    command(["minimap2", "-x", "sr", "--secondary=no", "-t", threads, gene, *reads],
            paf, prefix.with_suffix(".paf.log"))
    paf_stats(paf, name, length)  # Empty successful PAF is a valid zero-hit result.
    with paf.open("rb") as src, gzip.open(str(paf) + ".gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    paf.unlink()  # Only the newly generated uncompressed copy.
    bam = prefix.with_suffix(".bam")
    with prefix.with_suffix(".bam.log").open("wb") as log:
        mapping = subprocess.Popen(["minimap2", "-ax", "sr", "--secondary=no", "-t", str(threads), str(gene), *map(str, reads)],
                                   stdout=subprocess.PIPE, stderr=log)
        try:
            sorting = subprocess.run(["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"],
                                     stdin=mapping.stdout, stderr=log)
        finally:
            mapping.stdout.close()
        rc = mapping.wait()
        if rc or sorting.returncode:
            raise RuntimeError(f"Read mapping/sort failed: {prefix}; minimap2={rc}, sort={sorting.returncode}")
    subprocess.run(["samtools", "quickcheck", str(bam)], check=True)
    depth = prefix.with_suffix(".depth.tsv")
    command(["samtools", "depth", "-aa", bam], depth, prefix.with_suffix(".depth.log"))
    values = []
    with depth.open() as handle:
        for position, line in enumerate(handle, 1):
            ref, pos, value = line.rstrip().split("\t")
            if ref != name or int(pos) != position or int(value) < 0:
                raise ValueError("Invalid per-base depth/reference positions: " + str(depth))
            values.append(int(value))
    if len(values) != length:
        raise ValueError("Depth must explicitly include every reference base, including zeros: " + str(depth))
    bam.unlink()  # Only this helper's BAM; preserve FASTQs and original pipeline BAMs.
    return 100 * sum(v > 0 for v in values) / length


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wd", required=True, type=Path)
    parser.add_argument("--gene", required=True, type=input_file)
    parser.add_argument("--out", required=True, type=Path, help="Fixed managed output directory; replaced only after all selected rows succeed")
    parser.add_argument("--rows", default="all", help="Plan run_index IDs, e.g. 1-2 or 409-480,225-230")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sites", type=input_file)
    parser.add_argument("--catalogue", type=input_file)
    args = parser.parse_args()
    try:
        check_output(args.out)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    if bool(args.sites) != bool(args.catalogue):
        parser.error("Use --sites and --catalogue together")
    if args.threads < 1:
        parser.error("--threads must be positive")
    for tool in ("minimap2", "samtools"):
        if not shutil.which(tool):
            parser.error(f"{tool} unavailable in the active Linux environment")
    wd, gene = args.wd.resolve(), args.gene.resolve()
    header, sequence = fasta(gene)
    ref_name, length = header.split()[0], len(sequence)
    plan_path = wd / "input/rrmg/rrmg_permutation_plan.tsv"
    plan = table(plan_path)
    selected = None
    if args.rows != "all":
        selected = set()
        for part in args.rows.split(","):
            if not re.fullmatch(r"\d+(-\d+)?", part):
                parser.error("Invalid --rows")
            ends = [int(x) for x in part.split("-")]
            if ends[0] > ends[-1]:
                parser.error("Descending --rows range")
            selected.update(range(ends[0], ends[-1] + 1))
    indices = [int(r["run_index"]) for r in plan]
    if len(indices) != len(set(indices)) or len({r['run_id'] for r in plan}) != len(plan):
        parser.error("Duplicate plan IDs")
    if selected is not None and selected - set(indices):
        parser.error("Some --rows IDs are absent from plan")
    plan = [r for r in plan if selected is None or int(r["run_index"]) in selected]
    if not plan:
        parser.error("No selected plan rows")
    sites, site_audit = validated_sites(gene, args.sites, args.catalogue) if args.sites else ([], [])
    if site_audit:
        print(f"[SITES] {len(sites)} included; {sum(r['status'] == 'excluded' for r in site_audit)} excluded; see site_annotation_audit.tsv", flush=True)
    inputs = []
    # Check every sample before starting any mapping; do not silently drop failures.
    for row in plan:
        rid = row["run_id"]
        if Path(rid).name != rid or "\\" in rid or rid in (".", ".."):
            parser.error("Unsafe run_id")
        run = wd / "temp/rrmg_permutation" / rid
        contigs = run / "asm/contigs.fasta"
        if not contigs.is_file() or not contigs.stat().st_size:
            parser.error("Contigs unavailable: " + str(contigs))
        inputs.append((row, contigs, paired_reads(run / "qc", "*_1.fastq.gz"),
                       paired_reads(run / "reads", "*target_1.fastq.gz")))
    # All analysis inputs must live outside the replaceable output directory.
    source_files = [gene, plan_path, *([args.sites, args.catalogue] if args.sites else [])]
    source_files += [p for _, contigs, qc, target_reads in inputs for p in (contigs, *qc, *target_reads)]
    if any(p.resolve().is_relative_to(args.out.resolve()) for p in source_files):
        parser.error("--out must not contain analysis inputs")
    with managed_output(args.out) as out:
        target = out / "target_gene.fasta"
        shutil.copy2(gene, target)
        shutil.copy2(plan_path, out / "source_plan.tsv")
        if sites:
            write_table(out / "CYP51A_FungAMR_sites.tsv", sites, list(sites[0]))
            write_table(out / "site_annotation_audit.tsv", site_audit, ["label", "status", "reason"])
            shutil.copy2(args.sites, out / "source_sites.tsv")
            shutil.copy2(args.catalogue, out / "source_catalogue.tsv")
        manifest = []
        provenance = dict(gene_source=str(gene), reference_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                          catalogue_source=str(args.catalogue), sites_source=str(args.sites), runs=[])
        for row, contigs, qc, target_reads in inputs:
            plot_id = "row_" + str(int(row["run_index"])).zfill(5)
            output = out / plot_id
            output.mkdir()
            print("[MAP]", row["run_id"], flush=True)
            paf = output / "contigs_to_GENE.paf"
            command(["minimap2", "-x", "asm5", "--secondary=no", "-c", "-t", args.threads, target, contigs], paf, output / "contigs.log")
            _, contig_breadth = paf_stats(paf, ref_name, length)
            read_breadth = map_reads(target, qc, output / "reads_to_GENE", args.threads, ref_name, length)
            target_breadth = map_reads(target, target_reads, output / "target_only_reads_to_GENE", args.threads, ref_name, length)
            manifest.append(dict(PlotID=plot_id, RunID=row["run_id"], Seed=row["seed"],
                                 TargetAbundanceDesign=row["target_abundance"], TotalPairs=row["total_pairs"],
                                 TargetGenomeCoverageExpected=row.get("target_expected_cov", "NA") or "NA",
                                 TargetFA="target_gene.fasta", ContigPAF=f"{plot_id}/contigs_to_GENE.paf",
                                 ReadDepthTSV=f"{plot_id}/reads_to_GENE.depth.tsv", ReadPAFGZ=f"{plot_id}/reads_to_GENE.paf.gz",
                                 TargetOnlyReadDepthTSV=f"{plot_id}/target_only_reads_to_GENE.depth.tsv",
                                 TargetOnlyReadPAFGZ=f"{plot_id}/target_only_reads_to_GENE.paf.gz",
                                 ContigTargetBreadthPct=contig_breadth, ReadBreadthPct=read_breadth,
                                 TargetOnlyReadBreadthPct=target_breadth))
            provenance["runs"].append(dict(RunID=row["run_id"], inputs=[dict(path=str(p), size=p.stat().st_size,
                                                     mtime_ns=p.stat().st_mtime_ns) for p in (contigs, *qc, *target_reads)]))
        write_table(out / "plot_manifest.with_target_only.tsv", manifest, list(manifest[0]))
        # Marker only written after all requested rows and reference checks succeed.
        (out / "COMPLETE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"[DONE] {len(manifest)} runs; {len(sites)} validated site labels -> {args.out}")


if __name__ == "__main__":
    main()
