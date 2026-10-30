#!/usr/bin/env python3
"""Collect existing Step3/Step4 outputs for portable alignment figures.

This collector never runs minimap2 or samtools. It validates and copies the
contig PAF from Step3 plus mixed-read and target-only-read PAF/depth outputs
from Step4. Site annotations are validated against the supplied catalogue and
reference. No Step1-5 result is changed.
"""
import argparse
from contextlib import contextmanager
import csv
import gzip
import hashlib
import json
import re
import shutil
import tempfile
import uuid
from pathlib import Path
from rrmg_gene_evidence import (validated_evidence, validate_base_table, contig_details,
                                marker_path, mix_generation, generation_status)


def input_file(value):
    """Reject unset shell variables before Path('') can become the current directory."""
    if not value.strip():
        raise argparse.ArgumentTypeError("Empty file path; set SITES_SOURCE/CATALOGUE_SOURCE (runbook section 0 or 13).")
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


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else Path(path).open(encoding="utf-8")


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def paf_stats(path, name, length, min_mapq=0):
    covered = set()
    count = 0
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            try:
                valid = len(fields) >= 12 and fields[5] == name and int(fields[6]) == length
            except ValueError:
                valid = False
            if not valid:
                raise ValueError(f"Invalid PAF/reference mismatch: {path}")
            try:
                start, end = int(fields[7]), int(fields[8])
            except ValueError as exc:
                raise ValueError(f"Invalid PAF coordinates: {path}") from exc
            if not 0 <= start < end <= length:
                raise ValueError(f"PAF coordinate out of range: {path}")
            if int(fields[11]) < min_mapq:
                continue
            covered.update(range(start, end))
            count += 1
    return count, 100 * len(covered) / length


def depth_stats(path, name, length):
    values = []
    with Path(path).open(encoding="utf-8") as handle:
        for position, line in enumerate(handle, 1):
            fields = line.rstrip().split("\t")
            if len(fields) != 3:
                raise ValueError("Invalid per-base depth row: " + str(path))
            ref, pos, value = fields
            try:
                pos, value = int(pos), int(value)
            except ValueError as exc:
                raise ValueError("Non-integer depth/position: " + str(path)) from exc
            if ref != name or pos != position or value < 0:
                raise ValueError("Invalid per-base depth/reference positions: " + str(path))
            values.append(value)
    if len(values) != length:
        raise ValueError("Depth must explicitly include every reference base, including zeros: " + str(path))
    return 100 * sum(v > 0 for v in values) / length


def select_rows(plan, row_spec):
    selected = None
    if row_spec != "all":
        selected = set()
        for part in row_spec.split(","):
            if not re.fullmatch(r"\d+(-\d+)?", part):
                raise ValueError("Invalid --rows")
            ends = [int(x) for x in part.split("-")]
            if ends[0] > ends[-1]:
                raise ValueError("Descending --rows range")
            selected.update(range(ends[0], ends[-1] + 1))
    indices = [int(row["run_index"]) for row in plan]
    if len(indices) != len(set(indices)) or len({row["run_id"] for row in plan}) != len(plan):
        raise ValueError("Duplicate plan IDs")
    if selected is not None and selected - set(indices):
        raise ValueError("Some --rows IDs are absent from plan")
    result = [row for row in plan if selected is None or int(row["run_index"]) in selected]
    if not result:
        raise ValueError("No selected plan rows")
    return result


def require_input(path, allow_empty=False):
    path = Path(path)
    if not path.is_file() or (not allow_empty and path.stat().st_size == 0):
        raise ValueError("Missing/empty Step3/Step4 output: " + str(path))
    return path


def build_alignment_bundle(wd, gene, out, row_spec="all", sites_file=None, catalogue=None, plan_path=None,
                           require_base_evidence=False, root_path=None):
    """Build one alignment bundle from existing mappings; never run an aligner."""
    wd, gene, out = Path(wd).resolve(), Path(gene).resolve(), Path(out)
    plan_path = Path(plan_path).resolve() if plan_path else wd / "input/rrmg/rrmg_permutation_plan.tsv"
    header, sequence = fasta(gene)
    ref_name, length = header.split()[0], len(sequence)
    plan = select_rows(table(plan_path), row_spec)
    sites, site_audit = validated_sites(gene, sites_file, catalogue) if sites_file else ([], [])
    if site_audit:
        print(f"[SITES] {len(sites)} included; {sum(r['status'] == 'excluded' for r in site_audit)} excluded; see site_annotation_audit.tsv", flush=True)
    inputs = []
    # Check every sample before writing the bundle; do not silently drop failures.
    for row in plan:
        rid = row["run_id"]
        if Path(rid).name != rid or "\\" in rid or rid in (".", ".."):
            raise ValueError("Unsafe run_id")
        run = (Path(root_path) if root_path else wd / "temp/rrmg_permutation") / rid
        files = dict(
            contig_paf=require_input(run / "asm/contigs.fasta.to_GENE.paf", allow_empty=True),
            read_depth=require_input(run / "readmap/reads_to_GENE.depth.tsv"),
            read_paf=require_input(run / "readmap/reads_to_GENE.paf.gz"),
            target_depth=require_input(run / "readmap/target_only_reads_to_GENE.depth.tsv"),
            target_paf=require_input(run / "readmap/target_only_reads_to_GENE.paf.gz"),
        )
        _, contig_breadth = paf_stats(files["contig_paf"], ref_name, length, min_mapq=20)
        paf_stats(files["read_paf"], ref_name, length)
        paf_stats(files["target_paf"], ref_name, length)
        read_breadth = depth_stats(files["read_depth"], ref_name, length)
        target_breadth = depth_stats(files["target_depth"], ref_name, length)
        inputs.append((row, files, contig_breadth, read_breadth, target_breadth))

    source_files = [gene, plan_path, *([sites_file, catalogue] if sites_file else [])]
    source_files += [path for _, files, *_ in inputs for path in files.values()]
    if any(Path(path).resolve().is_relative_to(out.resolve()) for path in source_files):
        raise ValueError("--out must not contain analysis inputs")
    out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(gene, out / "target_gene.fasta")
    shutil.copy2(plan_path, out / "source_plan.tsv")
    if sites:
        write_table(out / "CYP51A_FungAMR_sites.tsv", sites, list(sites[0]))
        write_table(out / "site_annotation_audit.tsv", site_audit, ["label", "status", "reason"])
        shutil.copy2(sites_file, out / "source_sites.tsv")
        shutil.copy2(catalogue, out / "source_catalogue.tsv")
    manifest = []
    evidence_status, recoveries, differences = [], [], []
    reference_sha = hashlib.sha256(sequence.encode()).hexdigest()
    provenance = dict(
        collection_mode="reuse_existing_step3_step4_outputs",
        gene_source=str(gene),
        reference_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
        catalogue_source=str(catalogue or ""),
        sites_source=str(sites_file or ""),
        runs=[],
    )
    for row, files, contig_breadth, read_breadth, target_breadth in inputs:
        plot_id = "row_" + str(int(row["run_index"])).zfill(5)
        output = out / plot_id
        output.mkdir()
        print("[COLLECT]", row["run_id"], flush=True)
        destinations = dict(
            contig_paf=output / "contigs_to_GENE.paf",
            read_depth=output / "reads_to_GENE.depth.tsv",
            read_paf=output / "reads_to_GENE.paf.gz",
            target_depth=output / "target_only_reads_to_GENE.depth.tsv",
            target_paf=output / "target_only_reads_to_GENE.paf.gz",
        )
        for key, source in files.items():
            if key == "contig_paf":
                with source.open(encoding="utf-8") as inp, destinations[key].open("w", encoding="utf-8", newline="") as dst:
                    for line in inp:
                        fields = line.rstrip().split("\t")
                        if len(fields) >= 12 and int(fields[11]) >= 20:
                            dst.write(line)
            else:
                shutil.copy2(source, destinations[key])
        run = (Path(root_path) if root_path else wd / "temp/rrmg_permutation") / row["run_id"]
        current = mix_generation(run)
        evidence_paths = {"ContigCSPAF": run / "asm/contigs.fasta.to_GENE.cs.paf"}
        for key, prefix in (("Metagenomic", "reads_to_GENE"), ("Synthetic", "target_only_reads_to_GENE")):
            for level in ("raw", "q20"):
                evidence_paths[f"{key}_{level}"] = run / f"readmap/{prefix}.base_counts.{level}.tsv.gz"
        evidence = dict(RunID=row["run_id"])
        for key, path in evidence_paths.items():
            signature = None
            try:
                signature = validated_evidence(path, reference_sha, ref_name)
                if key != "ContigCSPAF":
                    validate_base_table(path, ref_name, sequence, signature)
                else:
                    paf_stats(path, ref_name, length)
            except (OSError, ValueError, KeyError) as exc:
                if require_base_evidence:
                    raise ValueError(f"Missing/invalid base evidence for {row['run_id']}: {path}; run rrmg_backfill_cyp51a_evidence.sh first ({exc})") from exc
                evidence[key] = "MISSING_OR_INVALID"
                continue
            evidence[key] = "PRESENT"
            shutil.copy2(path, output / path.name)
            shutil.copy2(marker_path(path), output / marker_path(path).name)
            generation = generation_status(current, signature.get("generation", ""))
            evidence[f"{key}_Generation"] = generation
            if key == "ContigCSPAF":
                recovery, changes = contig_details(path, sequence)
                recovery.update(RunID=row["run_id"], GenerationStatus=generation,
                                DiscordanceFlag="Y" if float(recovery["CombinedContigBreadth"]) >= 90 and
                                (read_breadth < 50 or float(recovery["CombinedContigBreadth"]) - read_breadth >= 50) else "N",
                                Interpretation="REVIEW_REQUIRED" if generation != "MATCH" else "ALIGNMENT_EVIDENCE_ONLY")
                recoveries.append(recovery)
                differences.extend(dict(RunID=row["run_id"], GenerationStatus=generation, **change) for change in changes)
        evidence_status.append(evidence)
        manifest.append(dict(
            PlotID=plot_id, RunID=row["run_id"], Seed=row["seed"],
            TargetAbundanceDesign=row["target_abundance"], TotalPairs=row["total_pairs"],
            TargetGenomeCoverageExpected=row.get("target_expected_cov", row.get("target_cov", "NA")) or "NA",
            TargetFA="target_gene.fasta", ContigPAF=f"{plot_id}/contigs_to_GENE.paf",
            ReadDepthTSV=f"{plot_id}/reads_to_GENE.depth.tsv", ReadPAFGZ=f"{plot_id}/reads_to_GENE.paf.gz",
            TargetOnlyReadDepthTSV=f"{plot_id}/target_only_reads_to_GENE.depth.tsv",
            TargetOnlyReadPAFGZ=f"{plot_id}/target_only_reads_to_GENE.paf.gz",
            ContigTargetBreadthPct=f"{contig_breadth:.6f}", ReadBreadthPct=f"{read_breadth:.6f}",
            TargetOnlyReadBreadthPct=f"{target_breadth:.6f}",
        ))
        provenance["runs"].append(dict(
            RunID=row["run_id"],
            inputs=[dict(path=str(path), size=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns,
                         sha256=digest(path)) for path in files.values()],
        ))
    write_table(out / "plot_manifest.with_target_only.tsv", manifest, list(manifest[0]))
    status_fields = ["RunID"]
    for key in evidence_paths:
        status_fields.extend((key, f"{key}_Generation"))
    write_table(out / "base_evidence_status.tsv", evidence_status, status_fields)
    write_table(out / "contig_recovery.tsv", recoveries,
                ["RunID", "RecoveryClass", "QualifyingContigN", "MaxSingleContigBreadth", "CombinedContigBreadth",
                 "DifferenceEventN", "GenerationStatus", "DiscordanceFlag", "Interpretation"])
    write_table(out / "contig_differences.tsv", differences,
                ["RunID", "ContigID", "GenePosition", "EventType", "RefAllele", "QueryAllele", "Strand", "MapQ", "Identity", "PassRecoveryFilter", "GenerationStatus"])
    provenance["base_evidence_complete"] = all(all(row.get(key) == "PRESENT" for key in evidence_paths) for row in evidence_status)
    provenance["recovery_filters"] = dict(min_mapq=20, min_block_length=150, min_identity_pct=95, full_breadth_pct=90)
    (out / "COMPLETE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return len(manifest), len(sites)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wd", required=True, type=Path)
    parser.add_argument("--gene", required=True, type=input_file)
    parser.add_argument("--out", required=True, type=Path, help="Fixed managed output directory")
    parser.add_argument("--rows", default="all", help="Plan run_index IDs, e.g. 1-2 or 409-480,225-230")
    parser.add_argument("--sites", type=input_file)
    parser.add_argument("--catalogue", type=input_file)
    parser.add_argument("--require-base-evidence", action="store_true", help="fail unless cs PAF and four base-count tables are valid for every run")
    args = parser.parse_args()
    try:
        check_output(args.out)
        if bool(args.sites) != bool(args.catalogue):
            raise ValueError("Use --sites and --catalogue together")
        with managed_output(args.out) as stage:
            # managed_output creates the staging directory; build into a child and
            # move its contents up so the ownership marker remains at --out.
            child = stage / "bundle"
            count, site_count = build_alignment_bundle(
                args.wd, args.gene, child, args.rows, args.sites, args.catalogue,
                require_base_evidence=args.require_base_evidence
            )
            for item in list(child.iterdir()):
                item.rename(stage / item.name)
            child.rmdir()
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f"[DONE] {count} runs; {site_count} validated site labels -> {args.out}")


if __name__ == "__main__":
    main()
