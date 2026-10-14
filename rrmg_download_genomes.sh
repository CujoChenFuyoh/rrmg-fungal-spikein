#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Download genomes for a flat RandomReadsMG community manifest.
#
# Input manifest columns:
#   genome_ID  name  taxid  rel_abund  fasta  role
#
# Compared with the earlier CAMISIM downloader, this does NOT write/update
# metadata.tsv, genome_to_id.tsv, or abundance.tsv. It only updates one simple
# RandomReadsMG manifest.
###############################################################################

usage() {
  cat <<'EOF'
Usage:
  rrmg_download_genomes.sh \
    --manifest input/rrmg/community_manifest.tsv \
    --out-manifest input/rrmg/community_manifest.downloaded.tsv \
    --threads 4

Required:
  --manifest FILE          Flat manifest: genome_ID, name, taxid, rel_abund, fasta, role

Optional:
  --out-manifest FILE      Output manifest after removing failed downloads and renormalizing abundance.
                           Default: overwrite input manifest in place.
  --genome-dir DIR         Override FASTA output directory. Default: dirname(manifest)/genomes
  --tmp-dir DIR            Temporary datasets directory. Default: dirname(manifest)/tmp_ncbi_datasets
  --success FILE           Success list. Default: dirname(manifest)/download_success.tsv
  --failed FILE            Failed list. Default: dirname(manifest)/download_failed.tsv
  --threads INT            Reserved for compatibility; downloads are currently sequential. Default: 4
  --keep-tmp               Do not remove per-genome temporary workdirs
  -h, --help               Show help

Dependencies:
  NCBI datasets CLI, unzip, python3
EOF
}

MANIFEST=""
OUT_MANIFEST=""
GENOME_DIR=""
TMP_DIR=""
SUCCESS_LIST=""
FAIL_LIST=""
THREADS=4
KEEP_TMP=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --out-manifest) OUT_MANIFEST="$2"; shift 2;;
    --genome-dir) GENOME_DIR="$2"; shift 2;;
    --tmp-dir) TMP_DIR="$2"; shift 2;;
    --success) SUCCESS_LIST="$2"; shift 2;;
    --failed) FAIL_LIST="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --keep-tmp) KEEP_TMP=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

: "${MANIFEST:?Missing --manifest}"
[[ -s "$MANIFEST" ]] || { echo "[ERROR] Manifest not found/empty: $MANIFEST" >&2; exit 1; }

BASE_DIR="$(cd "$(dirname "$MANIFEST")" && pwd)"
MANIFEST_ABS="$(cd "$(dirname "$MANIFEST")" && pwd)/$(basename "$MANIFEST")"
OUT_MANIFEST="${OUT_MANIFEST:-$MANIFEST_ABS}"
GENOME_DIR="${GENOME_DIR:-${BASE_DIR}/genomes}"
TMP_DIR="${TMP_DIR:-${BASE_DIR}/tmp_ncbi_datasets}"
SUCCESS_LIST="${SUCCESS_LIST:-${BASE_DIR}/download_success.tsv}"
FAIL_LIST="${FAIL_LIST:-${BASE_DIR}/download_failed.tsv}"
NORMALIZED_INPUT="${TMP_DIR}/manifest.normalized.tsv"

mkdir -p "$GENOME_DIR" "$TMP_DIR" "$BASE_DIR/download_logs"

command -v datasets >/dev/null 2>&1 || { echo "[ERROR] 'datasets' not found in PATH" >&2; exit 1; }
command -v unzip >/dev/null 2>&1 || { echo "[ERROR] 'unzip' not found in PATH" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "[ERROR] 'python3' not found in PATH" >&2; exit 1; }

cat <<EOF
[INFO] MANIFEST     : $MANIFEST_ABS
[INFO] OUT_MANIFEST : $OUT_MANIFEST
[INFO] GENOME_DIR   : $GENOME_DIR
[INFO] TMP_DIR      : $TMP_DIR
[INFO] THREADS      : $THREADS
EOF

# Normalize arbitrary manifest column spellings into a stable six-column TSV.
python3 - "$MANIFEST_ABS" "$NORMALIZED_INPUT" "$GENOME_DIR" <<'PY'
import csv, re, sys
from pathlib import Path

manifest = Path(sys.argv[1])
out = Path(sys.argv[2])
genome_dir = Path(sys.argv[3])

def norm(s):
    return re.sub(r"[\s_\-]+", "", str(s)).lower()

def pick(row, names, default=""):
    m = {norm(k): k for k in row}
    for n in names:
        k = norm(n)
        if k in m:
            return row.get(m[k], default)
    return default

def safe_id(s):
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s).strip())
    return s.strip("_") or "Genome"

with manifest.open(newline="") as fh:
    sample = fh.read(8192)
    fh.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
    except csv.Error:
        dialect = csv.excel_tab
    rows = list(csv.DictReader(fh, dialect=dialect))

if not rows:
    raise SystemExit("[ERROR] Empty manifest")

with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome_ID", "name", "taxid", "rel_abund", "fasta", "role"])
    for i, r in enumerate(rows, start=1):
        gid = pick(r, ["genome_ID", "genome", "id"]) or f"Genome{i}"
        gid = safe_id(gid)
        name = pick(r, ["name", "Name", "species", "taxon"]) or gid
        taxid = pick(r, ["taxid", "NCBI ID", "NCBI_ID", "TaxID"]).strip() or "NA"
        rel = pick(r, ["rel_abund", "relative_abundance", "abundance", "rel"], "0").strip() or "0"
        fasta = pick(r, ["fasta", "path_to_fasta_file", "genome_path", "path"])
        role = (pick(r, ["role"], "background") or "background").strip().lower()
        if not fasta:
            fasta = str(genome_dir / f"{gid}.fna")
        w.writerow([gid, name, taxid, rel, fasta, role])
PY

: > "$SUCCESS_LIST"
: > "$FAIL_LIST"

download_one_package () {
  local taxid="$1"
  local zipfile="$2"
  local mode="$3"      # reference / level
  local level="${4:-}" # complete/chromosome/scaffold/contig

  rm -f "$zipfile"

  if [[ "$mode" == "reference" ]]; then
    datasets download genome taxon "$taxid" \
      --reference \
      --include genome \
      --filename "$zipfile" \
      >"$BASE_DIR/download_logs/$(basename "$zipfile").log" 2>&1
  else
    datasets download genome taxon "$taxid" \
      --assembly-level "$level" \
      --include genome \
      --filename "$zipfile" \
      >"$BASE_DIR/download_logs/$(basename "$zipfile").log" 2>&1
  fi
}

pick_best_fna_from_package () {
  local workdir="$1"

  python3 - "$workdir" <<'PY'
import json, sys
from pathlib import Path

workdir = Path(sys.argv[1])
data_dir = workdir / "ncbi_dataset" / "data"
report = data_dir / "assembly_data_report.jsonl"

if not report.exists():
    print("")
    sys.exit(0)

records = []
with open(report, "r", encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        try:
            obj = json.loads(line)
        except Exception:
            continue

        accession = (
            obj.get("accession")
            or obj.get("assemblyInfo", {}).get("assemblyAccession")
            or obj.get("assembly", {}).get("accession")
            or ""
        )
        assembly_level = (
            obj.get("assembly_info", {}).get("assembly_level")
            or obj.get("assemblyInfo", {}).get("assemblyLevel")
            or obj.get("assemblyLevel")
            or ""
        )
        level_rank_map = {
            "Complete Genome": 4,
            "Chromosome": 3,
            "Scaffold": 2,
            "Contig": 1,
            "complete": 4,
            "chromosome": 3,
            "scaffold": 2,
            "contig": 1,
        }
        level_rank = level_rank_map.get(str(assembly_level), 0)
        refseq_rank = 1 if str(accession).startswith("GCF_") else 0
        records.append((refseq_rank, level_rank, accession))

if not records:
    print("")
    sys.exit(0)

records.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)

for _, _, acc in records:
    if not acc:
        continue
    subdir = data_dir / acc
    if not subdir.exists():
        continue
    candidates = sorted(list(subdir.glob("*_genomic.fna")) + list(subdir.glob("*.fna")) + list(subdir.glob("*.fa")) + list(subdir.glob("*.fasta")))
    for f in candidates:
        if f.is_file() and f.stat().st_size > 0:
            print(str(f))
            sys.exit(0)

print("")
PY
}

# Use the normalized manifest. Keep rows with existing FASTA without downloading.
tail -n +2 "$NORMALIZED_INPUT" | while IFS=$'\t' read -r genome_id name taxid rel_abund fasta role; do
  if [[ -z "$genome_id" ]]; then
    continue
  fi

  if [[ -s "$fasta" ]]; then
    echo "[SKIP] FASTA exists: $genome_id -> $fasta"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$genome_id" "$name" "$taxid" "$rel_abund" "$fasta" "$role" >> "$SUCCESS_LIST"
    continue
  fi

  if [[ -z "$taxid" || "$taxid" == "NA" || "$taxid" == "0" ]]; then
    echo "[FAIL] Missing usable taxid for $genome_id | $name"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$genome_id" "$name" "$taxid" "$rel_abund" "$fasta" "$role" >> "$FAIL_LIST"
    continue
  fi

  echo "[CHECK] $genome_id | taxid=$taxid | $name"

  workdir="${TMP_DIR}/${genome_id}"
  rm -rf "$workdir"
  mkdir -p "$workdir"

  success=0
  chosen_fna=""

  zipfile="${TMP_DIR}/${genome_id}.reference.zip"
  if download_one_package "$taxid" "$zipfile" "reference"; then
    unzip -q "$zipfile" -d "$workdir" || true
    datasets rehydrate --directory "$workdir" >/dev/null 2>&1 || true
    chosen_fna="$(pick_best_fna_from_package "$workdir")"
    if [[ -n "$chosen_fna" && -s "$chosen_fna" ]]; then
      success=1
      echo "[OK] reference genome selected for $genome_id"
    fi
  fi

  if [[ "$success" -eq 0 ]]; then
    for level in complete chromosome scaffold contig; do
      rm -rf "$workdir"
      mkdir -p "$workdir"
      zipfile="${TMP_DIR}/${genome_id}.${level}.zip"
      if download_one_package "$taxid" "$zipfile" "level" "$level"; then
        unzip -q "$zipfile" -d "$workdir" || true
        datasets rehydrate --directory "$workdir" >/dev/null 2>&1 || true
        chosen_fna="$(pick_best_fna_from_package "$workdir")"
        if [[ -n "$chosen_fna" && -s "$chosen_fna" ]]; then
          success=1
          echo "[OK] selected $level-level best available genome for $genome_id"
          break
        fi
      fi
    done
  fi

  if [[ "$success" -eq 1 && -n "$chosen_fna" && -s "$chosen_fna" ]]; then
    mkdir -p "$(dirname "$fasta")"
    cp -f "$chosen_fna" "$fasta"
    if [[ -s "$fasta" ]]; then
      printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$genome_id" "$name" "$taxid" "$rel_abund" "$fasta" "$role" >> "$SUCCESS_LIST"
      echo "[SAVE] $genome_id -> $fasta"
    else
      printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$genome_id" "$name" "$taxid" "$rel_abund" "$fasta" "$role" >> "$FAIL_LIST"
      echo "[FAIL] copy failed for $genome_id"
    fi
  else
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$genome_id" "$name" "$taxid" "$rel_abund" "$fasta" "$role" >> "$FAIL_LIST"
    echo "[FAIL] no usable genome found for $genome_id"
  fi

  if [[ "$KEEP_TMP" -eq 0 ]]; then
    rm -rf "$workdir"
  fi
done

python3 - "$SUCCESS_LIST" "$OUT_MANIFEST" "$FAIL_LIST" <<'PY'
import csv, sys
from pathlib import Path

success = Path(sys.argv[1])
out = Path(sys.argv[2])
rows = []
with open(sys.argv[3], newline="") as failed:
    if any(len(r) >= 6 and r[5].strip().lower() == "target"
           for r in csv.reader(failed, delimiter="\t")):
        raise SystemExit("[ERROR] Target genome download failed; inspect download_failed.tsv and download_logs. Output manifest not updated.")
if success.exists() and success.stat().st_size > 0:
    with success.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            gid, name, taxid, rel, fasta, role = parts[:6]
            try:
                rel_f = float(rel)
            except Exception:
                rel_f = 0.0
            rows.append({"genome_ID": gid, "name": name, "taxid": taxid, "rel_abund": rel_f, "fasta": fasta, "role": role or "background"})

if not rows:
    raise SystemExit("[ERROR] No genomes retained after download/check.")

# Renormalize positive rows. Zero-abundance target rows are preserved for later insertion.
pos_total = sum(max(0.0, r["rel_abund"]) for r in rows)
if pos_total > 0:
    for r in rows:
        if r["rel_abund"] > 0:
            r["rel_abund"] = r["rel_abund"] / pos_total

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome_ID", "name", "taxid", "rel_abund", "fasta", "role"])
    for r in rows:
        w.writerow([r["genome_ID"], r["name"], r["taxid"], f"{r['rel_abund']:.10f}", r["fasta"], r["role"]])

print(f"[OK] Updated manifest: {out}")
print(f"[OK] Retained genomes: {len(rows)}")
print(f"[OK] Positive abundance sum: {sum(r['rel_abund'] for r in rows):.10f}")
PY

echo "[INFO] Done."
echo "[INFO] Success count: $(wc -l < "$SUCCESS_LIST" 2>/dev/null || echo 0)"
if [[ -s "$FAIL_LIST" ]]; then
  echo "[INFO] Failed genomes saved to: $FAIL_LIST"
else
  echo "[INFO] No failed genomes."
fi
