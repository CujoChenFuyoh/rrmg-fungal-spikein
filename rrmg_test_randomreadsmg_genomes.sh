#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  rrmg_test_randomreadsmg_genomes.sh \
    --manifest input/rrmg/community_manifest.downloaded.tsv \
    --out-manifest input/rrmg/community_manifest.tested.tsv \
    -w /path/to/project \
    --test-pairs 200 --readlen 150 --seed 999

Purpose:
  Quick gatekeeper for FASTA compatibility with randomreadsmg.sh.
  Genomes that cannot generate a tiny test FASTQ are excluded and background
  abundances are renormalized. Target genomes are never silently removed.

Outputs:
  --out-manifest plus pass/fail/summary TSVs under --test-dir. Temporary test
  FASTQs are deleted after a successful check; per-genome logs are retained.

Required environment:
  RANDOMREADSMG_SH      Absolute path to randomreadsmg.sh
USAGE
}

MANIFEST=""
OUT_MANIFEST=""
WD="."
TEST_PAIRS=200
READLEN=150
SEED=999
TEST_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --out-manifest) OUT_MANIFEST="$2"; shift 2;;
    -w) WD="$2"; shift 2;;
    --test-pairs) TEST_PAIRS="$2"; shift 2;;
    --readlen) READLEN="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --test-dir) TEST_DIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

: "${MANIFEST:?Missing --manifest}"
: "${OUT_MANIFEST:?Missing --out-manifest}"
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/rrmg_environment.sh"
: "${RANDOMREADSMG_SH:?randomreadsmg.sh missing from PATH; set RANDOMREADSMG_SH}"

WD="$(mkdir -p "$WD" && cd "$WD" && pwd)"
cd "$WD"
[[ -s "$MANIFEST" ]] || { echo "[ERROR] Manifest not found: $MANIFEST" >&2; exit 1; }
[[ -s "$RANDOMREADSMG_SH" ]] || { echo "[ERROR] RANDOMREADSMG_SH not found: $RANDOMREADSMG_SH" >&2; exit 1; }

if [[ -z "$TEST_DIR" ]]; then
  TEST_DIR="$(dirname "$OUT_MANIFEST")/randomreadsmg_genome_test"
fi
mkdir -p "$TEST_DIR" "$TEST_DIR/logs" "$TEST_DIR/reads" "$(dirname "$OUT_MANIFEST")"

SUMMARY="$TEST_DIR/randomreadsmg_test_summary.tsv"
PASSED="$TEST_DIR/randomreadsmg_test_passed.tsv"
FAILED="$TEST_DIR/randomreadsmg_test_failed.tsv"

python3 - "$MANIFEST" "$OUT_MANIFEST" "$SUMMARY" "$PASSED" "$FAILED" "$TEST_DIR" "$RANDOMREADSMG_SH" "$TEST_PAIRS" "$READLEN" "$SEED" "$WD" <<'PY'
import csv
import os
import subprocess
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
out_manifest = Path(sys.argv[2])
summary = Path(sys.argv[3])
passed = Path(sys.argv[4])
failed = Path(sys.argv[5])
test_dir = Path(sys.argv[6])
rrmg = sys.argv[7]
test_pairs = int(float(sys.argv[8]))
readlen = int(float(sys.argv[9]))
base_seed = int(float(sys.argv[10]))
wd = Path(sys.argv[11])

def resolve(p):
    q = Path(p)
    return q if q.is_absolute() else wd / q

with manifest.open(newline="") as fh:
    sample = fh.read(8192); fh.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
    except csv.Error:
        dialect = csv.excel_tab
    rows = list(csv.DictReader(fh, dialect=dialect))
fields = rows[0].keys()

passed_rows = []
failed_rows = []
summary_rows = []

for i, row in enumerate(rows, 1):
    gid = row.get("genome_ID") or row.get("genome") or "Genome{}".format(i)
    fasta = row.get("fasta") or row.get("path_to_fasta_file") or row.get("path") or ""
    role = (row.get("role") or "background").lower()
    fa = resolve(fasta)
    out1 = test_dir / "reads" / "{}.test_1.fastq.gz".format(gid)
    out2 = test_dir / "reads" / "{}.test_2.fastq.gz".format(gid)
    log = test_dir / "logs" / "{}.randomreadsmg.test.log".format(gid)
    status = "FAIL"
    reason = ""
    if not fa.exists():
        reason = "FASTA_NOT_FOUND"
    else:
        cmd = ["bash", rrmg, "-Xmx" + os.environ.get("BBMAP_XMX", "2g"), "t=" + os.environ.get("THREADS", "1"), "overwrite=t", "ref={}".format(str(fa)), "out={}".format(str(out1)), "out2={}".format(str(out2)), "reads={}".format(test_pairs), "paired=t", "minlength={}".format(readlen), "maxlength={}".format(readlen), "seed={}".format(base_seed + i)]
        with log.open("w") as lfh:
            ret = subprocess.call(cmd, stdout=lfh, stderr=lfh)
        if ret == 0 and out1.exists() and out1.stat().st_size > 0 and out2.exists() and out2.stat().st_size > 0:
            status = "PASS"
            reason = "OK"
            try:
                out1.unlink(); out2.unlink()
            except Exception:
                pass
        else:
            reason = "RANDOMREADSMG_FAILED"
    srow = dict(row)
    srow.update({"test_status": status, "test_reason": reason, "test_log": str(log)})
    summary_rows.append(srow)
    if status == "PASS":
        passed_rows.append(row)
        print("[PASS] {}".format(gid), flush=True)
    else:
        failed_rows.append(srow)
        print("[FAIL] {} {}".format(gid, reason), flush=True)
        if role == "target":
            raise SystemExit("[ERROR] Target genome failed randomreadsmg test: {}".format(gid))

# renormalize background rows only; keep target rows unchanged
bg = [r for r in passed_rows if (r.get("role") or "background").lower() != "target"]
targets = [r for r in passed_rows if (r.get("role") or "background").lower() == "target"]
rel_field = "rel_abund" if "rel_abund" in fields else None
if rel_field:
    total = sum(float(r.get(rel_field) or 0) for r in bg)
    if total > 0:
        for r in bg:
            r[rel_field] = "{:.10f}".format(float(r.get(rel_field) or 0) / total)

for path, rws in [(out_manifest, bg + targets), (passed, passed_rows), (failed, failed_rows), (summary, summary_rows)]:
    if not rws:
        path.write_text("")
        continue
    fn = list(rws[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fn, delimiter="\t")
        w.writeheader(); w.writerows(rws)

print("[OK] wrote tested manifest: {}".format(out_manifest), flush=True)
print("[OK] passed: {} failed: {}".format(len(passed_rows), len(failed_rows)), flush=True)
PY
