#!/usr/bin/env python3
"""Pipeline stage 1B: collapse seed replicates into reusable background jobs.

The permutation plan contains one row per seed, but the background community
is intentionally identical across seeds. This script therefore creates one
background job per unique abundance/depth design and allocates integer read
pairs to background genomes with the largest-remainder method.

Inputs: permutation plan, referenced manifests, and all background FASTA files.
Outputs: ``background_plan.tsv``, ``background_index.tsv``, a FASTA-length
cache, and one ``background_design.tsv`` per unique background job.
"""

import csv
import gzip
import hashlib
import math
import re
import sys
from pathlib import Path


def usage():
    print("""Usage:
  python3 rrmg_make_background_plan.py \
    --plan input/rrmg/rrmg_permutation_plan.tsv \
    -w /path/to/project \
    --outdir temp/rrmg_background \
    --background-seed 1000

Purpose:
  Generate one background job for each unique design combination in rrmg_permutation_plan.tsv.
  Background is shared across seeds; only target reads vary by seed in Step 2.

Unique background key:
  manifest × target_abundance × total_pairs/depth × readlen × abundance_mode

Outputs:
  <outdir>/genome_length_cache.tsv
  <outdir>/background_plan.tsv
  <outdir>/background_index.tsv
  <outdir>/<bg_key>/background_design.tsv

Options:
  --length-cache FILE       Custom FASTA length cache path.
  --refresh-length-cache    Recalculate FASTA lengths.
  --min-pairs N             Drop background genomes assigned fewer than N read pairs. Default: 1.
  --readlen N               Default read length if plan lacks readlen. Default: 150.
""", file=sys.stderr)


def parse_args(argv):
    args = {
        "plan": "",
        "wd": "",
        "outdir": "temp/rrmg_background",
        "background_seed": 1000,
        "readlen": 150,
        "min_pairs": 1,
        "length_cache": "",
        "refresh_length_cache": False,
    }
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--plan":
            args["plan"] = argv[i + 1]; i += 2
        elif a == "-w":
            args["wd"] = argv[i + 1]; i += 2
        elif a == "--outdir":
            args["outdir"] = argv[i + 1]; i += 2
        elif a == "--background-seed":
            args["background_seed"] = int(float(argv[i + 1])); i += 2
        elif a == "--readlen":
            args["readlen"] = int(float(argv[i + 1])); i += 2
        elif a == "--min-pairs":
            args["min_pairs"] = int(float(argv[i + 1])); i += 2
        elif a == "--length-cache":
            args["length_cache"] = argv[i + 1]; i += 2
        elif a == "--refresh-length-cache":
            args["refresh_length_cache"] = True; i += 1
        elif a in ("-h", "--help"):
            usage(); sys.exit(0)
        else:
            raise SystemExit("[ERROR] Unknown arg: {}".format(a))
    if not args["plan"]:
        raise SystemExit("[ERROR] Missing --plan")
    if not args["wd"]:
        raise SystemExit("[ERROR] Missing -w")
    return args


def safe(s):
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s).strip())
    return s.strip("_") or "NA"


def norm(s):
    return re.sub(r"[\s_\-]+", "", str(s)).lower()


def resolve(wd, p):
    if not p:
        return ""
    q = Path(str(p))
    if not q.is_absolute():
        q = wd / q
    return str(q)


def open_text(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(str(path), "rt")
    return path.open("r", errors="ignore")


def fasta_len(path):
    n = 0
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n


def sniff_dialect(path):
    with Path(path).open("r", errors="ignore") as fh:
        sample = fh.read(8192)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;")
    except csv.Error:
        return csv.excel_tab


def read_manifest(wd, manifest):
    manifest = Path(resolve(wd, manifest))
    if not manifest.exists():
        raise SystemExit("[ERROR] Manifest not found: {}".format(manifest))
    dialect = sniff_dialect(manifest)
    with manifest.open(newline="", errors="ignore") as fh:
        rows = list(csv.DictReader(fh, dialect=dialect))
    out = []
    for i, r in enumerate(rows, 1):
        keys = {norm(k): k for k in r.keys()}
        def get(*names):
            for name in names:
                k = norm(name)
                if k in keys:
                    return r.get(keys[k], "")
            return ""
        gid = get("genome_ID", "genome", "id") or "Genome{}".format(i)
        name = get("name", "species", "taxon") or gid
        taxid = get("taxid", "NCBI ID", "NCBI_ID", "TaxID") or "NA"
        fasta = resolve(wd, get("fasta", "path_to_fasta_file", "genome_path", "path"))
        role = (get("role") or "background").lower()
        try:
            rel = float(get("rel_abund", "relative_abundance", "abundance", "rel") or 0)
        except Exception:
            rel = 0.0
        out.append({"gid": gid, "name": name, "taxid": taxid, "fasta": fasta, "role": role, "rel": max(0.0, rel)})
    return out


def is_target(row, target_name, target_taxid, target_gid):
    if row.get("role", "").lower() == "target":
        return True
    if target_gid and row.get("gid") == target_gid:
        return True
    if target_taxid and target_taxid not in {"NA", "0", ""} and str(row.get("taxid")) == str(target_taxid):
        return True
    if target_name and row.get("name", "").strip().lower() == target_name.strip().lower():
        return True
    return False


def read_plan(path):
    with Path(path).open(newline="", errors="ignore") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit("[ERROR] Empty plan: {}".format(path))
    return rows


def read_length_cache(path):
    cache = {}
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return cache
    with path.open(newline="", errors="ignore") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            fasta = r.get("fasta", "")
            length = r.get("length_bp", "")
            if not fasta or not length:
                continue
            try:
                cache[fasta] = int(float(length))
            except Exception:
                pass
    return cache


def write_length_cache(path, cache):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fasta", "length_bp"], delimiter="\t")
        w.writeheader()
        for fasta in sorted(cache):
            w.writerow({"fasta": fasta, "length_bp": cache[fasta]})
    tmp.replace(path)


def largest_remainder(weights, total):
    total = int(round(float(total)))
    if total <= 0 or not weights:
        return [0] * len(weights)
    s = sum(weights)
    if s <= 0:
        return [0] * len(weights)
    raw = [w / s * total for w in weights]
    base = [int(math.floor(x)) for x in raw]
    rem = total - sum(base)
    order = sorted(range(len(raw)), key=lambda i: raw[i] - base[i], reverse=True)
    for i in order[:rem]:
        base[i] += 1
    return base


def get_float(row, key, default=0.0):
    v = row.get(key, "")
    if v is None or v == "":
        return default
    return float(v)


def collect_needed_fastas(wd, plan_rows):
    manifests = {}
    fastas = set()
    for row in plan_rows:
        manifest = resolve(wd, row.get("manifest", ""))
        if manifest not in manifests:
            manifests[manifest] = read_manifest(wd, manifest)
        target_name = row.get("target_name", "")
        target_taxid = row.get("target_taxid", "")
        target_gid = row.get("target_genome_id", "")
        for r in manifests[manifest]:
            if is_target(r, target_name, target_taxid, target_gid):
                continue
            if r["rel"] <= 0:
                continue
            if not r["fasta"]:
                raise SystemExit("[ERROR] Missing FASTA path for {}".format(r["gid"]))
            if not Path(r["fasta"]).exists():
                raise SystemExit("[ERROR] FASTA not found for {}: {}".format(r["gid"], r["fasta"]))
            fastas.add(r["fasta"])
    return manifests, sorted(fastas)


def ensure_lengths(fastas, cache, refresh=False):
    total = len(fastas)
    updated = False
    for i, fasta in enumerate(fastas, 1):
        if (not refresh) and fasta in cache and cache[fasta] > 0:
            print("[LEN] {}/{} cached {}".format(i, total, fasta), flush=True)
            continue
        print("[LEN] {}/{} scanning {}".format(i, total, fasta), flush=True)
        length = fasta_len(fasta)
        if length <= 0:
            raise SystemExit("[ERROR] FASTA has zero length: {}".format(fasta))
        cache[fasta] = length
        updated = True
    return updated


def total_pairs_from_row(row, readlen, bg_rows):
    total_pairs = row.get("total_pairs", "")
    total_bases = row.get("total_bases", "")
    total_cov = row.get("total_cov", "")
    community_size = row.get("community_size_bp", "")
    if total_pairs:
        return int(round(float(total_pairs)))
    if total_bases:
        return int(round(float(total_bases) / (readlen * 2)))
    if total_cov:
        if community_size:
            csize = float(community_size)
        else:
            csize = sum(r["bg_weight"] * r["length"] for r in bg_rows)
        return int(round(float(total_cov) * csize / (readlen * 2)))
    raise SystemExit("[ERROR] Plan row has no total_pairs/total_bases/total_cov")


def background_key(row, manifest_abs, target_abundance, total_pairs_i, readlen, abundance_mode):
    sample = row.get("sample") or row.get("run_id") or "bg"
    digest = hashlib.md5("{}|{}|{}|{}|{}".format(manifest_abs, target_abundance, total_pairs_i, readlen, abundance_mode).encode("utf-8")).hexdigest()[:8]
    return "{}__bg_{}".format(safe(sample), digest)


def main():
    args = parse_args(sys.argv[1:])
    wd = Path(args["wd"]).resolve()
    plan_path = Path(resolve(wd, args["plan"]))
    outdir = Path(resolve(wd, args["outdir"]))
    outdir.mkdir(parents=True, exist_ok=True)
    length_cache_path = Path(resolve(wd, args["length_cache"])) if args["length_cache"] else outdir / "genome_length_cache.tsv"
    bg_plan = outdir / "background_plan.tsv"
    bg_index = outdir / "background_index.tsv"

    print("[INFO] wd={}".format(wd), flush=True)
    print("[INFO] plan={}".format(plan_path), flush=True)
    print("[INFO] outdir={}".format(outdir), flush=True)
    print("[INFO] length_cache={}".format(length_cache_path), flush=True)

    plan_rows = read_plan(plan_path)
    manifests, fastas = collect_needed_fastas(wd, plan_rows)
    print("[INFO] unique manifests: {}".format(len(manifests)), flush=True)
    print("[INFO] unique background FASTAs: {}".format(len(fastas)), flush=True)
    cache = read_length_cache(length_cache_path)
    print("[INFO] cached FASTA lengths: {}".format(len(cache)), flush=True)
    updated = ensure_lengths(fastas, cache, refresh=args["refresh_length_cache"])
    if updated or args["refresh_length_cache"] or (not length_cache_path.exists()):
        write_length_cache(length_cache_path, cache)
        print("[OK] wrote length cache: {}".format(length_cache_path), flush=True)
    else:
        print("[OK] using existing length cache", flush=True)

    jobs = {}
    for row_i, row in enumerate(plan_rows, 1):
        if row_i == 1 or row_i == len(plan_rows) or row_i % 10 == 0:
            print("[PLAN] processing row {}/{}".format(row_i, len(plan_rows)), flush=True)
        manifest_abs = resolve(wd, row.get("manifest", ""))
        target_abundance = get_float(row, "target_abundance", 0.0)
        readlen = int(float(row.get("readlen") or args["readlen"]))
        abundance_mode = row.get("abundance_mode", "read_fraction") or "read_fraction"
        target_name = row.get("target_name", "")
        target_taxid = row.get("target_taxid", "")
        target_gid = row.get("target_genome_id", "")
        mrows = manifests[manifest_abs]
        bg = [dict(r) for r in mrows if (not is_target(r, target_name, target_taxid, target_gid)) and r["rel"] > 0]
        if not bg:
            raise SystemExit("[ERROR] No background genomes left for target {}".format(target_name))
        bg_total = sum(r["rel"] for r in bg)
        if bg_total <= 0:
            raise SystemExit("[ERROR] Background total relative abundance is zero")
        for r in bg:
            r["bg_weight"] = r["rel"] / bg_total
            r["length"] = cache[r["fasta"]]
        total_pairs_i = total_pairs_from_row(row, readlen, bg)
        target_pairs = int(round(total_pairs_i * target_abundance))
        bg_pairs_total = max(0, total_pairs_i - target_pairs)
        bg_key = background_key(row, manifest_abs, target_abundance, total_pairs_i, readlen, abundance_mode)
        if bg_key in jobs:
            continue
        pairs = largest_remainder([r["bg_weight"] for r in bg], bg_pairs_total)
        bg_dir = outdir / bg_key
        reads_dir = bg_dir / "reads"
        design = bg_dir / "background_design.tsv"
        bg_r1 = reads_dir / "{}.background_1.fastq.gz".format(bg_key)
        bg_r2 = reads_dir / "{}.background_2.fastq.gz".format(bg_key)
        reads_dir.mkdir(parents=True, exist_ok=True)
        with design.open("w", newline="") as dfh:
            w = csv.writer(dfh, delimiter="\t")
            w.writerow(["idx", "genome_ID", "name", "taxid", "role", "rel_abund_final", "bg_weight", "fasta", "genome_len_bp", "sim_pairs", "sim_cov", "seed"])
            idx = 0
            for r, p in zip(bg, pairs):
                if p < args["min_pairs"]:
                    continue
                idx += 1
                seed = args["background_seed"] + idx
                sim_cov = p * readlen * 2 / r["length"]
                rel_final = r["bg_weight"] * (1.0 - target_abundance)
                w.writerow([idx, r["gid"], r["name"], r["taxid"], "background", "{:.10f}".format(rel_final), "{:.10f}".format(r["bg_weight"]), r["fasta"], r["length"], p, "{:.10f}".format(sim_cov), seed])
        jobs[bg_key] = {
            "bg_index": len(jobs) + 1,
            "bg_key": bg_key,
            "manifest": manifest_abs,
            "target_abundance": "{:.10f}".format(target_abundance),
            "total_pairs": str(total_pairs_i),
            "readlen": str(readlen),
            "background_pairs": str(bg_pairs_total),
            "bg_dir": str(bg_dir),
            "design_tsv": str(design),
            "bg_r1": str(bg_r1),
            "bg_r2": str(bg_r2),
            "status": "pending",
        }
    fields = ["bg_index", "bg_key", "manifest", "target_abundance", "total_pairs", "readlen", "background_pairs", "bg_dir", "design_tsv", "bg_r1", "bg_r2", "status"]
    for path in (bg_plan, bg_index):
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for k in sorted(jobs, key=lambda x: jobs[x]["bg_index"]):
                w.writerow(jobs[k])
    print("[OK] Background jobs: {}".format(len(jobs)), flush=True)
    print("[OK] Wrote: {}".format(bg_plan), flush=True)
    print("[OK] Wrote: {}".format(bg_index), flush=True)
    print("[OK] Length cache: {}".format(length_cache_path), flush=True)
    print("[OK] Submit Step 1 array as: qsub -J 1-{} rrmg_step1_background_array.pbs".format(len(jobs)), flush=True)


if __name__ == "__main__":
    main()
