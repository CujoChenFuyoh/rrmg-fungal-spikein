#!/usr/bin/env bash
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_ROOT="$(mktemp -d)"
trap 'rm -rf -- "$TEST_ROOT"' EXIT

PIPELINE="$TEST_ROOT/pipeline"
WD="$TEST_ROOT/work"
FAKEBIN="$TEST_ROOT/bin"
JAVA_HOME_FAKE="$TEST_ROOT/java"
CALLS="$TEST_ROOT/calls"
CONFIG="$TEST_ROOT/rrmg_config.sh"
QSUB_LOG="$TEST_ROOT/qsub.log"

mkdir -p \
  "$PIPELINE" \
  "$WD/input/rrmg" \
  "$WD/temp/rrmg_background" \
  "$WD/ref" \
  "$FAKEBIN" \
  "$JAVA_HOME_FAKE/bin" \
  "$CALLS"

for helper in rrmg_run_one.sh rrmg_environment.sh; do
  cp "$REPO/$helper" "$PIPELINE/$helper"
done

printf '%s\n' \
  $'run_index\trun_id' \
  $'1\trun_1' \
  $'2\trun_2' \
  $'3\trun_3' \
  $'4\trun_4' \
  > "$WD/input/rrmg/rrmg_permutation_plan.tsv"

printf '%s\n' \
  $'bg_index\tbg_key' \
  $'1\tbg_1' \
  $'2\tbg_2' \
  $'3\tbg_3' \
  > "$WD/temp/rrmg_background/background_plan.tsv"

printf '%s\n' $'bg_key\tbg_r1\tbg_r2' $'bg_1\ta\tb' \
  > "$WD/temp/rrmg_background/background_index.tsv"
printf '>gene\nACGT\n' > "$WD/ref/CYP51A_wt.fasta"
printf '>target\nACGTACGT\n' > "$WD/ref/Afum_Af293.fasta"

printf '%s\n' \
  '#!/usr/bin/env bash' \
  'conda() {' \
  '  [[ "${1:-}" == "activate" ]] && return 0' \
  '  return 0' \
  '}' \
  > "$TEST_ROOT/conda.sh"

printf '%s\n' \
  '#!/usr/bin/env bash' \
  'exit 0' \
  > "$TEST_ROOT/fake_command.sh"
chmod +x "$TEST_ROOT/fake_command.sh"

for command_name in minimap2 samtools vamb randomreadsmg.sh randomreads.sh; do
  cp "$TEST_ROOT/fake_command.sh" "$FAKEBIN/$command_name"
done
cp "$TEST_ROOT/fake_command.sh" "$JAVA_HOME_FAKE/bin/java"
chmod +x "$FAKEBIN"/* "$JAVA_HOME_FAKE/bin/java"

printf '%s\n' \
  '#!/usr/bin/env bash' \
  'set -euo pipefail' \
  'name="$(basename "$0")"' \
  'row="${RRMG_ROW_INDEX:?missing RRMG_ROW_INDEX}"' \
  '{' \
  '  printf "step=%s\nrow=%s\n" "$name" "$row"' \
  '  printf "threads=%s\nbbmap_xmx=%s\n" "${THREADS:-}" "${BBMAP_XMX:-}"' \
  '  printf "args="' \
  '  printf "%q " "$@"' \
  '  printf "\n"' \
  '} > "${SMOKE_OUT:?}/${name}.${row}.txt"' \
  'if [[ "${SMOKE_FAIL_STEP:-}" == "$name" && "${SMOKE_FAIL_ROW:-}" == "$row" ]]; then' \
  '  exit 9' \
  'fi' \
  > "$TEST_ROOT/fake_step.sh"
chmod +x "$TEST_ROOT/fake_step.sh"

STEP_SCRIPTS=(
  rrmg_step1_make_background_reads.sh
  rrmg_step2_make_target_and_merge.sh
  rrmg_step3_downstream_analysis.sh
  rrmg_step3b_vamb_binning.sh
  rrmg_step4_5_supplementary_downstream.sh
  rrmg_step4_reads_to_gene_depth.sh
  rrmg_step5_afum_target_reads_genome_cov.sh
)
for step_script in "${STEP_SCRIPTS[@]}"; do
  cp "$TEST_ROOT/fake_step.sh" "$PIPELINE/$step_script"
done
cp "$TEST_ROOT/fake_step.sh" "$PIPELINE/rrmg_evaluate_vamb_bins.py"

{
  printf 'PIPELINE_DIR=%q\n' "$PIPELINE"
  printf 'WD=%q\n' "$WD"
  printf 'CONDA_SH=%q\n' "$TEST_ROOT/conda.sh"
  printf 'RRMG_CONDA_ENV=%q\n' 'fungamr_meta_smoke'
  printf 'BBMAP_DIR=%q\n' "$FAKEBIN"
  printf 'RANDOMREADSMG_SH=%q\n' "$FAKEBIN/randomreadsmg.sh"
  printf 'RANDOMREADS_SH=%q\n' "$FAKEBIN/randomreads.sh"
  printf 'JAVA_HOME=%q\n' "$JAVA_HOME_FAKE"
  printf 'JAVA_BIN=%q\n' "$JAVA_HOME_FAKE/bin/java"
  printf 'PLAN=%q\n' 'input/rrmg/rrmg_permutation_plan.tsv'
  printf 'BACKGROUND_PLAN=%q\n' 'temp/rrmg_background/background_plan.tsv'
  printf 'BACKGROUND_INDEX=%q\n' 'temp/rrmg_background/background_index.tsv'
  printf 'GENE_REF=%q\n' 'ref/CYP51A_wt.fasta'
  printf 'TARGET_GENOME=%q\n' 'ref/Afum_Af293.fasta'
  printf 'RUN_KRAKEN=0\n'
  printf 'VAMB_COMMAND=vamb\n'
  printf 'STEP2_KEEP_TARGET_READS=1\n'
  printf 'STEP45_MAKE_BINS=1\n'
  printf 'LOCAL_JOBS=1\n'
} > "$CONFIG"

run_local() {
  PATH="$FAKEBIN:$PATH" \
  RRMG_CONFIG="$CONFIG" \
  SMOKE_OUT="$CALLS" \
  SMOKE_FAIL_STEP="${SMOKE_FAIL_STEP:-}" \
  SMOKE_FAIL_ROW="${SMOKE_FAIL_ROW:-}" \
    bash "$REPO/${SMOKE_RUNNER:-rrmg_submit.sh}" "$@"
}

assert_file() {
  [[ -s "$1" ]] || { echo "[TEST ERROR] Expected file missing: $1" >&2; exit 1; }
}

run_local step1 1-2 --threads 3 --memory 7 --jobs 2
assert_file "$CALLS/rrmg_step1_make_background_reads.sh.1.txt"
assert_file "$CALLS/rrmg_step1_make_background_reads.sh.2.txt"
grep -q '^threads=3$' "$CALLS/rrmg_step1_make_background_reads.sh.1.txt"
grep -q '^bbmap_xmx=7g$' "$CALLS/rrmg_step1_make_background_reads.sh.1.txt"

run_local step2 2-4 -t 1 -m 6 --jobs 2
for row in 2 3 4; do
  assert_file "$CALLS/rrmg_step2_make_target_and_merge.sh.${row}.txt"
done
grep -q -- '--keep-target-reads' "$CALLS/rrmg_step2_make_target_and_merge.sh.2.txt"
grep -q '^threads=1$' "$CALLS/rrmg_step2_make_target_and_merge.sh.2.txt"
grep -q '^bbmap_xmx=6g$' "$CALLS/rrmg_step2_make_target_and_merge.sh.2.txt"

run_local step3 1 -t 5 -m 37
run_local step3b 1 -t 6
run_local step45 1 -t 7
run_local step4 2 -t 8
run_local step5 3 -t 9
grep -q -- "--skip-target-genome" "$CALLS/rrmg_step4_5_supplementary_downstream.sh.2.txt"
grep -q -- "--skip-read-to-gene" "$CALLS/rrmg_step4_5_supplementary_downstream.sh.3.txt"
assert_file "$CALLS/rrmg_step3_downstream_analysis.sh.1.txt"
assert_file "$CALLS/rrmg_step3b_vamb_binning.sh.1.txt"
assert_file "$CALLS/rrmg_step4_5_supplementary_downstream.sh.1.txt"
grep -q -- '-t 5 -mem 37' "$CALLS/rrmg_step3_downstream_analysis.sh.1.txt"
grep -q -- '--only-assembler' "$CALLS/rrmg_step3_downstream_analysis.sh.1.txt"

# Config controls the same native mode for all entry points, even old configs
# without the new key. Bad values fail instead of silently changing the method.
printf 'STEP3_ONLY_ASSEMBLER=0\n' >> "$CONFIG"
run_local step3 1
grep -q -- '--with-error-correction' "$CALLS/rrmg_step3_downstream_analysis.sh.1.txt"
if grep -q -- '--only-assembler' "$CALLS/rrmg_step3_downstream_analysis.sh.1.txt"; then exit 1; fi
printf 'STEP3_ONLY_ASSEMBLER=invalid\n' >> "$CONFIG"
if run_local step3 1 >/dev/null 2>&1; then
  echo '[TEST ERROR] Invalid assembly mode unexpectedly succeeded' >&2
  exit 1
fi
printf 'STEP3_ONLY_ASSEMBLER=1\n' >> "$CONFIG"
grep -q -- '-t 6' "$CALLS/rrmg_step3b_vamb_binning.sh.1.txt"
grep -q '^threads=7$' "$CALLS/rrmg_step4_5_supplementary_downstream.sh.1.txt"
grep -q '^threads=8$' "$CALLS/rrmg_step4_5_supplementary_downstream.sh.2.txt"
grep -q '^threads=9$' "$CALLS/rrmg_step4_5_supplementary_downstream.sh.3.txt"

run_local step2 1 --force
grep -q -- '--force' "$CALLS/rrmg_step2_make_target_and_merge.sh.1.txt"

if run_local step2 1-5 --dry-run >/dev/null 2>&1; then
  echo '[TEST ERROR] Out-of-range selection unexpectedly succeeded' >&2
  exit 1
fi
if run_local step3 1 --threads 0 >/dev/null 2>&1; then
  echo '[TEST ERROR] Zero threads unexpectedly succeeded' >&2
  exit 1
fi
if run_local step3 1 --memory 0 >/dev/null 2>&1; then
  echo '[TEST ERROR] Zero memory unexpectedly succeeded' >&2
  exit 1
fi
if run_local step5 1 --memory 8 >/dev/null 2>&1; then
  echo '[TEST ERROR] Unsupported Step5 memory override unexpectedly succeeded' >&2
  exit 1
fi

if SMOKE_FAIL_STEP=rrmg_step3_downstream_analysis.sh \
   SMOKE_FAIL_ROW=2 \
   run_local step3 1-2 --jobs 2 >/dev/null 2>&1; then
  echo '[TEST ERROR] Injected row failure was not propagated' >&2
  exit 1
fi

# Any accidental scheduler call must fail. No .pbs files exist in the pipeline.
printf '#!/bin/bash\nexit 99\n' > "$FAKEBIN/qsub"
chmod +x "$FAKEBIN/qsub"
run_local step2 2-3 --jobs 1
SMOKE_RUNNER=rrmg_run_local.sh run_local step2 1 --dry-run
if run_local step2 1 -W depend=afterok:123 >/dev/null 2>&1; then
  echo '[TEST ERROR] Old qsub option was accepted' >&2
  exit 1
fi
# A parent scheduler's array ID must not override an explicitly selected row.
PBS_ARRAY_INDEX=999 run_local step2 4
# A sparse original plan is selected by ID, not by physical line number.
printf 'run_index\trun_id\n409\trun409\n410\trun410\n' > "$WD/input/rrmg/rrmg_permutation_plan.tsv"
run_local step2 409-410
assert_file "$CALLS/rrmg_step2_make_target_and_merge.sh.410.txt"
if run_local step2 1 --dry-run >/dev/null 2>&1; then exit 1; fi
printf 'run_index\trun_id\n1\trun1\n1\trun2\n' > "$WD/input/rrmg/rrmg_permutation_plan.tsv"
if run_local step2 1 --dry-run >/dev/null 2>&1; then exit 1; fi
echo '[OK] Both scheduler-free entry points, all stages, force, concurrency, failure and ID selection.'
