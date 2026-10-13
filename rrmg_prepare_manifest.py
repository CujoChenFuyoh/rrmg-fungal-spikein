#!/usr/bin/env python3
"""
Pipeline stage 0A: prepare the community manifest used by all later steps.

Input
-----
One CSV/TSV taxon table (for example ``result/sample_taxa_table.tsv``), or an
existing manifest when ``--input-is-manifest`` is supplied.

Output
------
``community_manifest.tsv`` with the stable columns
``genome_ID, name, taxid, rel_abund, fasta, role``.

Prepare a simple RandomReadsMG / BBMap randomreads.sh community manifest.

This script replaces the CAMISIM-style metadata.tsv/genome_to_id.tsv/abundance.tsv
input stack with one flat TSV:

  genome_ID\tname\ttaxid\trel_abund\tfasta\trole

role is normally "background". A target species can be added as role="target"
with rel_abund=0 for downloading; the actual target abundance is usually set
later by rrmg_virtual_metagenome_permutation.sh or a permutation plan.

No pandas dependency.

Important: unless ``--target-abundance`` is explicitly supplied, the target
row is kept at abundance 0. Its abundance is injected later by
``rrmg_make_plan.py`` / Step 2, so it does not contaminate the background.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def sniff_dialect(path: Path) -> csv.Dialect:
    sample = path.read_text(errors="ignore")[:8192]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;")
    except csv.Error:
        class Tab(csv.Dialect):
            delimiter = "\t"
            quotechar = '"'
            doublequote = True
            skipinitialspace = False
            lineterminator = "\n"
            quoting = csv.QUOTE_MINIMAL
        return Tab


def norm_key(s: str) -> str:
    return re.sub(r"[\s_\-]+", "", str(s)).lower()


def pick(row: Dict[str, str], names: Iterable[str], default: str = "") -> str:
    for n in names:
        if n in row:
            return row.get(n, default)
    nmap = {norm_key(k): k for k in row.keys()}
    for n in names:
        k = norm_key(n)
        if k in nmap:
            return row.get(nmap[k], default)
    return default


def to_float(x: object, default: float = 0.0) -> float:
    try:
        s = str(x).strip()
        if s == "":
            return default
        return float(s)
    except Exception:
        return default


def to_taxid(x: object) -> str:
    s = str(x).strip()
    if not s:
        return "NA"
    try:
        return str(int(float(s)))
    except Exception:
        return s


def read_table(path: Path) -> List[Dict[str, str]]:
    dialect = sniff_dialect(path)
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, dialect=dialect)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"[ERROR] Empty input table: {path}")
    return rows


def safe_id(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s.strip())
    return s.strip("_") or "Genome"


def parse_rows(args: argparse.Namespace) -> List[Dict[str, object]]:
    raw_rows = read_table(Path(args.input))

    include_re: Optional[re.Pattern[str]] = re.compile(args.include_name_regex) if args.include_name_regex else None
    exclude_re: Optional[re.Pattern[str]] = re.compile(args.exclude_name_regex) if args.exclude_name_regex else None

    items: List[Dict[str, object]] = []
    input_is_manifest = args.input_is_manifest

    for row in raw_rows:
        if input_is_manifest:
            gid = pick(row, [args.genome_id_column, "genome_ID", "genome", "id"])
            name = pick(row, [args.name_column, "name", "Name", "species", "taxon"]).strip() or gid
            taxid = to_taxid(pick(row, [args.taxid_column, "taxid", "NCBI ID", "NCBI_ID", "TaxID"]))
            rank = pick(row, [args.rank_column, "rank", "NCBI Rank", "NCBI_Rank"]).strip()
            abundance = to_float(pick(row, [args.abundance_column, "rel_abund", "relative_abundance", "abundance", "rel"]), 0.0)
            fasta = pick(row, [args.fasta_column or "fasta", "path_to_fasta_file", "genome_path", "path"])
            role = pick(row, ["role"], "background") or "background"
        else:
            gid = ""
            name = pick(row, [args.name_column, "Name", "name", "taxon", "species"]).strip()
            taxid = to_taxid(pick(row, [args.taxid_column, "NCBI ID", "NCBI_ID", "taxid", "TaxID"]))
            rank = pick(row, [args.rank_column, "NCBI Rank", "NCBI_Rank", "rank"]).strip()
            abundance = to_float(pick(row, [args.abundance_column, "Summed base count", "Base count", "rel_abund", "relative_abundance", "abundance"]), 0.0)
            fasta = pick(row, [args.fasta_column], "") if args.fasta_column else ""
            role = "background"

        if rank and args.rank and rank.lower() != args.rank.lower():
            continue
        if taxid in {"", "0", "0.0", "NA", "nan", "None"}:
            continue
        if include_re and not include_re.search(name):
            continue
        if exclude_re and exclude_re.search(name):
            continue
        if abundance < 0:
            continue
        if not input_is_manifest and abundance <= 0:
            continue

        items.append({
            "genome_ID": gid,
            "name": name,
            "taxid": taxid,
            "abundance_raw": abundance,
            "fasta": fasta,
            "role": role,
        })

    items.sort(key=lambda x: float(x["abundance_raw"]), reverse=True)
    if args.topN and args.topN > 0:
        # Apply topN to background rows only. Target will be added after this.
        bg = [x for x in items if str(x.get("role", "background")).lower() != "target"]
        non_bg = [x for x in items if str(x.get("role", "background")).lower() == "target"]
        items = bg[: args.topN] + non_bg

    if not items:
        raise SystemExit("[ERROR] No rows left after filtering.")

    # Assign genome IDs and FASTA paths.
    used_ids = set()
    for i, item in enumerate(items, start=1):
        gid = str(item.get("genome_ID") or "").strip()
        if not gid:
            gid = f"{args.genome_prefix}{i}"
        base_gid = safe_id(gid)
        gid = base_gid
        suffix = 2
        while gid in used_ids:
            gid = f"{base_gid}_{suffix}"
            suffix += 1
        used_ids.add(gid)
        item["genome_ID"] = gid

        if not str(item.get("fasta") or "").strip():
            item["fasta"] = args.fasta_template.format(
                outdir=str(Path(args.outdir)),
                genome_ID=gid,
                name=safe_id(str(item["name"])),
                taxid=item["taxid"],
            )

    return items


def target_matches(item: Dict[str, object], target_name: str, target_taxid: str, target_gid: str) -> bool:
    if str(item.get("role", "")).lower() == "target":
        return True
    if target_gid and str(item.get("genome_ID", "")) == target_gid:
        return True
    if target_taxid and target_taxid not in {"NA", ""} and str(item.get("taxid", "")) == target_taxid:
        return True
    if target_name and str(item.get("name", "")).strip().lower() == target_name.strip().lower():
        return True
    return False


def add_or_mark_target(items: List[Dict[str, object]], args: argparse.Namespace) -> List[Dict[str, object]]:
    target_name = args.target_name or ""
    target_taxid = to_taxid(args.target_taxid) if args.target_taxid else ""
    target_gid = args.target_genome_id or ""
    target_fasta = args.target_fasta or ""

    if not (target_name or target_taxid or target_gid or target_fasta):
        return items

    matches = [x for x in items if target_matches(x, target_name, target_taxid, target_gid)]
    if matches:
        target = matches[0]
        target["role"] = "target"
        if target_name:
            target["name"] = target_name
        if target_taxid:
            target["taxid"] = target_taxid
        if target_gid:
            target["genome_ID"] = target_gid
        if target_fasta:
            target["fasta"] = target_fasta
        if args.target_abundance is None:
            # Keep the row available for downloading, but do not let it contaminate
            # the background community unless user explicitly sets abundance here.
            target["abundance_raw"] = 0.0
    else:
        if not target_name:
            target_name = target_gid or f"TaxID_{target_taxid}"
        if not target_taxid:
            target_taxid = "NA"
        if not target_gid:
            target_gid = safe_id(target_name)
        if not target_fasta:
            target_fasta = args.fasta_template.format(
                outdir=str(Path(args.outdir)),
                genome_ID=target_gid,
                name=safe_id(target_name),
                taxid=target_taxid,
            )
        items.append({
            "genome_ID": target_gid,
            "name": target_name,
            "taxid": target_taxid,
            "abundance_raw": 0.0 if args.target_abundance is None else float(args.target_abundance),
            "fasta": target_fasta,
            "role": "target",
        })

    return items


def normalize(items: List[Dict[str, object]], target_abundance: Optional[float]) -> List[Dict[str, object]]:
    target_rows = [x for x in items if str(x.get("role", "background")).lower() == "target"]
    background = [x for x in items if str(x.get("role", "background")).lower() != "target"]

    bg_total = sum(max(0.0, float(x.get("abundance_raw", 0.0))) for x in background)
    if bg_total <= 0:
        raise SystemExit("[ERROR] Background abundance total is 0 after filtering.")

    if target_abundance is None:
        target_abundance = None
        bg_scale = 1.0
    else:
        if not (0 <= target_abundance < 1):
            raise SystemExit(f"[ERROR] --target-abundance must be >=0 and <1, got {target_abundance}")
        bg_scale = 1.0 - target_abundance

    out: List[Dict[str, object]] = []
    for x in background:
        rel = max(0.0, float(x.get("abundance_raw", 0.0))) / bg_total * bg_scale
        y = dict(x)
        y["rel_abund"] = rel
        y["role"] = "background"
        out.append(y)

    if target_rows:
        t = dict(target_rows[0])
        if target_abundance is None:
            t["rel_abund"] = max(0.0, float(t.get("abundance_raw", 0.0)))
            # If target was added only for downloading, rel_abund stays 0.
        else:
            t["rel_abund"] = target_abundance
        t["role"] = "target"
        out.append(t)

    # Final correction for positive rows only. Preserve 0-abundance target download row.
    positive_total = sum(float(x.get("rel_abund", 0.0)) for x in out)
    if positive_total > 0:
        for x in out:
            if float(x.get("rel_abund", 0.0)) > 0:
                x["rel_abund"] = float(x["rel_abund"]) / positive_total

    return out


def write_manifest(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["genome_ID", "name", "taxid", "rel_abund", "fasta", "role"])
        for x in rows:
            writer.writerow([
                x["genome_ID"],
                x["name"],
                x["taxid"],
                f"{float(x.get('rel_abund', 0.0)):.10f}",
                x["fasta"],
                x.get("role", "background"),
            ])


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Prepare one flat RandomReadsMG community manifest, optionally with a target species row."
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Input CSV/TSV taxa table or existing manifest")
    ap.add_argument("-o", "--outdir", required=True,
                    help="Output directory")
    ap.add_argument("--output", default=None,
                    help="Output manifest path. Default: <outdir>/community_manifest.tsv")
    ap.add_argument("--input-is-manifest", action="store_true",
                    help="Input already has manifest-like columns")

    ap.add_argument("--topN", type=int, default=30,
                    help="Keep top N background taxa after filtering; use 0 to keep all. Default: 30")
    ap.add_argument("--rank", default="species",
                    help="Keep only this NCBI Rank if rank column exists. Use empty string to disable. Default: species")
    ap.add_argument("--exclude-name-regex", default=None)
    ap.add_argument("--include-name-regex", default=None)

    ap.add_argument("--abundance-column", default="Summed base count")
    ap.add_argument("--name-column", default="Name")
    ap.add_argument("--taxid-column", default="NCBI ID")
    ap.add_argument("--rank-column", default="NCBI Rank")
    ap.add_argument("--genome-id-column", default="genome_ID")
    ap.add_argument("--fasta-column", default=None,
                    help="Optional FASTA path column in input")
    ap.add_argument("--fasta-template", default="{outdir}/genomes/{genome_ID}.fna",
                    help="Template for FASTA paths. Fields: outdir, genome_ID, name, taxid")
    ap.add_argument("--genome-prefix", default="Genome")

    ap.add_argument("--target-name", default=None,
                    help="Target species to insert/mark, e.g. 'Aspergillus fumigatus'")
    ap.add_argument("--target-taxid", default=None,
                    help="Target NCBI TaxID, e.g. 746128 or 5085 depending on desired taxon")
    ap.add_argument("--target-genome-id", default=None,
                    help="Genome ID to use for target row. Default derived from target name")
    ap.add_argument("--target-fasta", default=None,
                    help="Target genome FASTA path. If omitted, uses fasta-template")
    ap.add_argument("--target-abundance", type=float, default=None,
                    help="Optional fixed target abundance for this manifest. Usually omit and set during permutation.")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "genomes").mkdir(parents=True, exist_ok=True)

    rows = parse_rows(args)
    rows = add_or_mark_target(rows, args)
    rows = normalize(rows, args.target_abundance)

    out_path = Path(args.output) if args.output else outdir / "community_manifest.tsv"
    write_manifest(out_path, rows)

    n_bg = sum(1 for x in rows if str(x.get("role", "background")).lower() != "target")
    n_target = sum(1 for x in rows if str(x.get("role", "background")).lower() == "target")
    total = sum(float(x.get("rel_abund", 0.0)) for x in rows)
    print(f"[OK] Wrote manifest: {out_path}")
    print(f"[OK] Background rows: {n_bg}; target rows: {n_target}; rel_abund sum: {total:.10f}")
    print("[OK] Columns: genome_ID, name, taxid, rel_abund, fasta, role")


if __name__ == "__main__":
    main()
