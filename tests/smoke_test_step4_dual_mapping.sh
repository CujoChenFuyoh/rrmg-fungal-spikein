#!/usr/bin/env bash
# Optional real-tool integration test for both Step4 read sets.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for tool in minimap2 samtools python3; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[SKIP] $tool is required for the Step4 dual-mapping integration test."
    exit 0
  fi
done

ROOT="$(mktemp -d "${TMPDIR:-/tmp}/rrmg-step4-dual.XXXXXX")"
trap 'rm -rf -- "$ROOT"' EXIT
WD="$ROOT/work"
RUN="$WD/temp/rrmg_permutation/Pilot__Seed1"
mkdir -p "$WD/input/rrmg" "$WD/ref" "$RUN/design" "$RUN/qc" "$RUN/reads"

python3 - "$WD" <<'PY'
import gzip
import json
import random
import sys
from pathlib import Path

wd = Path(sys.argv[1])
run = wd / "temp/rrmg_permutation/Pilot__Seed1"
rng = random.Random(1729)
gene = "".join(rng.choice("ACGT") for _ in range(650))
(wd / "ref/gene.fasta").write_text(">CYP51A_toy\n" + gene + "\n", encoding="utf-8")
(wd / "input/rrmg/plan.tsv").write_text(
    "run_index\trun_id\tsample\ttarget_abundance\ttotal_pairs\tseed\ttarget_genome_id\n"
    "1\tPilot__Seed1\tPilot\t0.4\t10\t1\tTarget\n",
    encoding="utf-8",
)

def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def pair(name, start):
    r1 = gene[start:start + 75]
    r2 = reverse_complement(gene[start + 125:start + 200])
    return (
        f"@{name}/1\n{r1}\n+\n{'I' * len(r1)}\n",
        f"@{name}/2\n{r2}\n+\n{'I' * len(r2)}\n",
    )

target = [pair(f"target{i}", start) for i, start in enumerate((10, 100, 200, 350), 1)]
background = []
for i in range(6):
    r1 = "".join(rng.choice("ACGT") for _ in range(75))
    r2 = "".join(rng.choice("ACGT") for _ in range(75))
    background.append((f"@background{i}/1\n{r1}\n+\n{'I' * 75}\n",
                       f"@background{i}/2\n{r2}\n+\n{'I' * 75}\n"))

def write_pair(prefix, records):
    for mate in (0, 1):
        path = Path(str(prefix) + f"_{mate + 1}.fastq.gz")
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as out:
                for record in records:
                    out.write(record[mate].encode())

write_pair(run / "reads/Target.Seed1.target", target)
write_pair(run / "qc/Pilot.Seed1", background + target)
generation = "0123456789abcdef0123456789abcdef"
(run / "design/read_mix_audit.tsv").write_text(
    "key\tvalue\n"
    f"mix_generation_id\t{generation}\n"
    "total_pairs_written\t10\n"
    "background_pairs_written\t6\n"
    "target_pairs_written\t4\n"
    "gzip_members_per_output\t1\n",
    encoding="utf-8",
)
(run / "qc/Pilot.Seed1.fastp.json").write_text(
    json.dumps({"summary": {"before_filtering": {"total_reads": 20}}}) + "\n",
    encoding="utf-8",
)
(run / "qc/Pilot.Seed1.mix_generation.txt").write_text(generation + "\n", encoding="utf-8")
PY

THREADS=2 bash "$REPO/rrmg_step4_reads_to_gene_depth.sh" \
  "$WD/input/rrmg/plan.tsv" 1 -w "$WD" -gene "$WD/ref/gene.fasta"

for prefix in reads_to_GENE target_only_reads_to_GENE; do
  for suffix in coverage_summary.tsv depth.tsv paf.gz; do
    test -s "$RUN/readmap/${prefix}.${suffix}"
  done
  gzip -t "$RUN/readmap/${prefix}.paf.gz"
done

python3 - "$RUN" <<'PY'
import csv
import gzip
import sys
from pathlib import Path

run = Path(sys.argv[1])
for prefix in ("reads_to_GENE", "target_only_reads_to_GENE"):
    with (run / f"readmap/{prefix}.coverage_summary.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1 and int(rows[0]["MappedReadAlignments"]) > 0, (prefix, rows)
    assert float(rows[0]["BreadthPct"]) > 0, (prefix, rows)
    with gzip.open(run / f"readmap/{prefix}.paf.gz", "rt") as handle:
        assert any(line.strip() for line in handle), prefix
PY

MIXED_HASH="$(sha256sum "$RUN/readmap/reads_to_GENE.coverage_summary.tsv" | awk '{print $1}')"
rm -f "$RUN/readmap/target_only_reads_to_GENE."*
THREADS=2 bash "$REPO/rrmg_step4_reads_to_gene_depth.sh" \
  "$WD/input/rrmg/plan.tsv" 1 -w "$WD" -gene "$WD/ref/gene.fasta"
test "$MIXED_HASH" = "$(sha256sum "$RUN/readmap/reads_to_GENE.coverage_summary.tsv" | awk '{print $1}')"
test -s "$RUN/readmap/target_only_reads_to_GENE.coverage_summary.tsv"

echo "[OK] Step4 mapped mixed and target-only reads and resumed each set independently."
