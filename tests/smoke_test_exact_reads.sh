#!/usr/bin/env bash
# Optional local integration test: real BBMap RandomReads -> exact mixer -> fastp.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RANDOMREADS_SH="${RANDOMREADS_SH:-${BBMAP_DIR:-}/randomreads.sh}"
if [[ ! -s "$RANDOMREADS_SH" ]] || ! command -v fastp >/dev/null 2>&1; then
  echo "[SKIP] Set RANDOMREADS_SH (or BBMAP_DIR) and put fastp on PATH for the integration test."
  exit 0
fi

ROOT="$(mktemp -d "${TMPDIR:-/tmp}/rrmg-exact-reads.XXXXXX")"
trap 'rm -rf -- "$ROOT"' EXIT
WD="$ROOT/work"
mkdir -p "$WD/input/rrmg" "$WD/ref" "$WD/temp/background/reads" "$WD/temp/rrmg_background"

python3 - "$WD" <<'PY'
import gzip
from pathlib import Path
import sys

root = Path(sys.argv[1])
with (root / "ref/target.fasta").open("w") as handle:
    for i in range(8):
        sequence = ("ACGTGCACTG"[(i % 10):] + "TGCACGTACG"[:(i % 10)]) * 300
        handle.write(f">chr{i + 1}\n{sequence}\n")

def record(name, mate):
    seq = ("ACGT" if mate == 1 else "TGCA") * 38
    seq = seq[:150]
    return f"@{name}/{mate}\n{seq}\n+\n{'I' * 150}\n".encode()

for mate in (1, 2):
    path = root / f"temp/background/reads/background_{mate}.fastq.gz"
    for start, stop in ((0, 40), (40, 100)):
        with path.open("ab") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as out:
                for i in range(start, stop):
                    out.write(record(f"background{i}", mate))
PY

BG1="$WD/temp/background/reads/background_1.fastq.gz"
BG2="$WD/temp/background/reads/background_2.fastq.gz"
printf 'key\tvalue\nsource_pairs\t100\ndesign_pairs\t90\nr1\t%s\nr2\t%s\n' "$BG1" "$BG2" > "${BG1}.pair_count.tsv"
printf 'genome_ID\tname\ttaxid\trel_abund\tfasta\trole\nTarget\tTarget species\t1\t0\t%s\ttarget\n' \
  "$WD/ref/target.fasta" > "$WD/input/rrmg/manifest.tsv"
printf 'run_index\trun_id\tsample\tmanifest\ttarget_abundance\ttotal_pairs\treadlen\tseed\ttarget_name\ttarget_taxid\ttarget_genome_id\ttarget_fasta\n1\tPilot__Seed1\tPilot\t%s\t0.3076923076923077\t130\t150\t1\tTarget species\t1\tTarget\t%s\n' \
  "$WD/input/rrmg/manifest.tsv" "$WD/ref/target.fasta" > "$WD/input/rrmg/plan.tsv"
printf 'bg_key\tmanifest\ttarget_abundance\ttotal_pairs\treadlen\tbg_r1\tbg_r2\nbackground\t%s\t0.3076923076923077\t130\t150\t%s\t%s\n' \
  "$WD/input/rrmg/manifest.tsv" "$BG1" "$BG2" > "$WD/temp/rrmg_background/background_index.tsv"

RANDOMREADS_SH="$RANDOMREADS_SH" THREADS=1 BBMAP_XMX=1g \
  bash "$REPO/rrmg_step2_make_target_and_merge.sh" \
  "$WD/input/rrmg/plan.tsv" 1 -w "$WD" \
  --background-index "$WD/temp/rrmg_background/background_index.tsv" --keep-target-reads

RUN="$WD/temp/rrmg_permutation/Pilot__Seed1"
TARGET1="$RUN/reads/Target.Seed1.target_1.fastq.gz"
TARGET2="$RUN/reads/Target.Seed1.target_2.fastq.gz"
MIX1="$RUN/reads/Pilot.Seed1_1.fastq.gz"
MIX2="$RUN/reads/Pilot.Seed1_2.fastq.gz"
AUDIT="$RUN/design/read_mix_audit.tsv"
[[ "$(python3 "$REPO/rrmg_mix_paired_fastq.py" count --r1 "$TARGET1" --r2 "$TARGET2")" == 40 ]]
python3 "$REPO/rrmg_mix_paired_fastq.py" validate-mix \
  --audit "$AUDIT" --expected-pairs 130 --r1 "$MIX1" --r2 "$MIX2"

mkdir -p "$RUN/qc" "$RUN/asm"
QC1="$RUN/qc/Pilot.Seed1_1.fastq.gz"
QC2="$RUN/qc/Pilot.Seed1_2.fastq.gz"
FASTP_JSON="$RUN/qc/Pilot.Seed1.fastp.json"
fastp -w 2 -i "$MIX1" -I "$MIX2" \
  -o "$QC1" -O "$QC2" -j "$FASTP_JSON" -h "$RUN/qc/Pilot.Seed1.fastp.html" \
  >/dev/null 2>"$RUN/qc/Pilot.Seed1.fastp.log"
python3 "$REPO/rrmg_mix_paired_fastq.py" validate-fastp --json "$FASTP_JSON" --expected-pairs 130
[[ "$(python3 "$REPO/rrmg_mix_paired_fastq.py" count --r1 "$MIX1" --r2 "$MIX2")" == 130 ]]
MIX_GENERATION="$(awk -F '\t' '$1=="mix_generation_id"{print $2}' "$AUDIT")"
printf '%s\n' "$MIX_GENERATION" > "$RUN/qc/Pilot.Seed1.mix_generation.txt"

sleep 1
cp "$WD/ref/target.fasta" "$RUN/asm/contigs.fasta"
printf '%s\n' "$MIX_GENERATION" > "$RUN/asm/RRMG_MIX_GENERATION"
printf '%s\n' 'only-assembler' > "$RUN/asm/RRMG_ASSEMBLY_MODE"
bash "$REPO/rrmg_step3_downstream_analysis.sh" "$WD/input/rrmg/plan.tsv" 1 \
  -w "$WD" -gene "$WD/ref/target.fasta" --skip-kraken --skip-assemble -t 2 -mem 1
THREADS=2 bash "$REPO/rrmg_step4_reads_to_gene_depth.sh" "$WD/input/rrmg/plan.tsv" 1 \
  -w "$WD" -gene "$WD/ref/target.fasta"
THREADS=2 bash "$REPO/rrmg_step5_afum_target_reads_genome_cov.sh" "$WD/input/rrmg/plan.tsv" 1 \
  -w "$WD" --target-genome "$WD/ref/target.fasta"

sleep 1
RANDOMREADS_SH="$RANDOMREADS_SH" THREADS=1 BBMAP_XMX=1g \
  bash "$REPO/rrmg_step2_make_target_and_merge.sh" \
  "$WD/input/rrmg/plan.tsv" 1 -w "$WD" \
  --background-index "$WD/temp/rrmg_background/background_index.tsv" --keep-target-reads --force
if bash "$REPO/rrmg_step3_downstream_analysis.sh" "$WD/input/rrmg/plan.tsv" 1 \
     -w "$WD" -gene "$WD/ref/target.fasta" --skip-kraken --skip-assemble -t 2 -mem 1 \
     >"$ROOT/stale_step3.log" 2>&1; then
  echo "[TEST ERROR] Step3 reused outputs after Step2 was rebuilt" >&2
  exit 1
fi
grep -q 'older Step 2 read-mix generation' "$ROOT/stale_step3.log"
if THREADS=2 bash "$REPO/rrmg_step4_reads_to_gene_depth.sh" "$WD/input/rrmg/plan.tsv" 1 \
     -w "$WD" -gene "$WD/ref/target.fasta" >"$ROOT/stale_step4.log" 2>&1; then
  echo "[TEST ERROR] Step4 reused QC after Step2 was rebuilt" >&2
  exit 1
fi
grep -q 'older/unaudited Step 2 mix' "$ROOT/stale_step4.log"

echo "[OK] Exact Step2, full fastp input, Step3/4/5 audits and stale-output refusal all passed."
