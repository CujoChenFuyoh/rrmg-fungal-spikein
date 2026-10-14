# RRMG fungal spike-in permutation workflow

Version: `0.3.0` · [GitHub](https://github.com/CujoChenFuyoh/rrmg-fungal-spikein) · [中文教程](docs/INSTALL.zh-CN.md)

Simulate a reusable metagenomic background, spike in a fungal target across
abundance × depth × seed combinations, then perform QC, assembly, gene mapping
and read coverage analysis. One-sample VAMB binning is optional.

> **Changed in 0.3.0:** `rrmg_submit.sh` now runs selected rows directly. It never
> calls `qsub` or executes PBS files. Activate your environment first. Scheduling
> belongs to your own job script; [examples](examples/README.md) cover local,
> ordinary PBS and PBS arrays. Old qsub options passed to this command are rejected.

> **Scientific limitation:** Step3 writes an empty VCF placeholder and `VarN=0`.
> It measures gene recovery, NOT validated fungAMR mutation calls. Collectors
> merge existing results; a successful collection does NOT prove every planned
> replicate has finished. Read-count semantics are unchanged by this release.

## Workflow and files

| Stage | Script | Required input | Output |
|---|---|---|---|
| Prepare manifest | `rrmg_prepare_manifest.py` | Taxon abundance TSV/CSV | Six-column manifest |
| Download genomes | `rrmg_download_genomes.sh` | Manifest with taxids | FASTA, downloaded manifest, success/failure lists |
| Test genomes | `rrmg_test_randomreadsmg_genomes.sh` | Downloaded manifest, FASTA, BBMap | Tested manifest and per-genome test logs |
| Plan | `rrmg_make_plan.py` | Tested manifest, target genome, design parameters | Permutation plan |
| Background plan | `rrmg_make_background_plan.py` | Permutation plan, genomes | Background plan, index and per-condition design |
| step1 | `rrmg_step1_make_background_reads.sh` | Background plan/design, genomes | Reusable background FASTQ pairs |
| step2 | `rrmg_step2_make_target_and_merge.sh` | Plan, background index/FASTQ, target genome | Target-only reads, mixed reads, merge design |
| step3 | `rrmg_step3_downstream_analysis.sh` | Mixed reads, gene FASTA | fastp QC, metaSPAdes contigs, contig-to-gene mapping |
| step4 | `rrmg_step4_reads_to_gene_depth.sh` | Plan/design, QC reads, gene FASTA | Gene depth and coverage summary |
| step5 | `rrmg_step5_afum_target_reads_genome_cov.sh` | Plan/design, target-only reads, target genome | Genome coverage and optional coverage windows |
| step45 | `rrmg_step4_5_supplementary_downstream.sh` | Inputs for step4 + step5 | Both supplementary analyses |
| step3b (optional) | `rrmg_step3b_vamb_binning.sh` | QC reads, contigs, target genome from design | VAMB bins and recovery metrics |
| Collect | `rrmg_collect_*.py`, `rrmg_collect_supplementary_downstream.sh` | Existing per-run results | Merged TSVs |

Step5 coverage windows are not metagenomic bins. VAMB is a separate optional step.
Step3 QC uses **fastp**, not FastQC. Enable Kraken2 only when its database is available.

```text
SpikeIN/
├── rrmg_config.sh                # private deployment config
├── input/rrmg/                  # manifests, genomes, permutation plan
├── ref/                         # target genome and target gene
├── temp/rrmg_background/        # shared background design, index, reads
├── temp/rrmg_permutation/<run_id>/
│   ├── design/
│   ├── reads/
│   ├── qc/
│   ├── asm/contigs.fasta
│   ├── map/
│   ├── readmap/
│   ├── target_readmap/
│   └── binning/                 # optional VAMB
├── result/
└── logs/local/<invocation>/     # row_*.log and summary.tsv (also inside PBS)
```

## 1. Install and activate

Linux/Ubuntu (including WSL), Bash and Python >=3.10 are required for analysis.
Use the existing integrated environment. `datasets`, Java and analysis tools
are all used from that environment; no separate ncbi_datasets environment is needed.

```bash
conda activate fungamr_meta
git clone https://github.com/CujoChenFuyoh/rrmg-fungal-spikein.git
cd rrmg-fungal-spikein
chmod +x rrmg_submit.sh rrmg_run_local.sh
python -m pip install .
rrmg --version
```

Or install without a checkout:

```bash
python -m pip install 'git+https://github.com/CujoChenFuyoh/rrmg-fungal-spikein.git@main'
rrmg init --workdir "$HOME/SpikeIN"
```

Without `--pipeline-dir`, init exports a versioned `.rrmg/core-0.3.0/` snapshot.
Existing configs/data are never overwritten. To pin code, replace `main` with
a tested, pushed commit SHA. Local edits are installed by `pip install .`, not
a Git URL. This package is not automatically published to PyPI.
`pip` installs code, not bioinformatics binaries or experimental data.
If creating a new environment, use `conda env create -f environment.yml` first.
Save `conda list --explicit` after validation for exact Linux environment reuse.

## 2. Configure a new local project

From the checkout, with `fungamr_meta` active:

```bash
export CODE="$PWD"
export WD="$HOME/SpikeIN_local_test"
rrmg init --workdir "$WD" --pipeline-dir "$CODE"
export RRMG_CONFIG="$WD/rrmg_config.sh"
nano "$RRMG_CONFIG"
mkdir -p "$WD/input/rrmg" "$WD/ref" "$WD/temp" "$WD/result"
rrmg doctor --config "$RRMG_CONFIG"
```

Only `PIPELINE_DIR` and `WD` must be site-specific. Tool overrides normally stay
empty: Java and `randomreadsmg.sh` resolve from the active PATH. If BBMap is
installed separately, set `RANDOMREADSMG_SH=/absolute/path/bbmap/randomreadsmg.sh`
in your private config (or `BBMAP_DIR` for PATH lookup). No environment change is needed.
Some PBS systems need explicit Java selection; see the commented block in both
[PBS examples](examples/README.md). Release configs contain no machine's Java path.

The default up-to-4 threads / 4g simulation heap / 16 GB SPAdes memory are pilot
settings, not a guarantee for a large metagenome. Adjust per-row settings to
your hardware/job allocation, leaving memory overhead. Step2 uses one thread
by default for seeded simulation reproducibility. `--jobs N` multiplies
these requirements. Missing future plan/reference files in doctor are warnings;
missing programs are errors. Doctor does not download inputs.

## 3. Prepare manifest → download genomes → compatibility test

**New experiments only.** For existing results, keep the original manifest, plan,
IDs, seeds, genomes and background index. Do not recreate them to resume.

Taxon input requires columns `Name`, `NCBI ID`, `NCBI Rank`, `Summed base count`.
Different column names can be supplied via `rrmg_prepare_manifest.py --help`.
The included `examples/taxa_demo.tsv` is a two-species demonstration, not research data.
For real analyses use your own abundance table. Default `--topN` is 30; here it
is explicitly 2 for the demonstration.

```bash
source "$RRMG_CONFIG"
python3 "$PIPELINE_DIR/rrmg_prepare_manifest.py" \
  --input "$PIPELINE_DIR/examples/taxa_demo.tsv" \
  --outdir "$WD/input/rrmg" --topN 2 \
  --target-name 'Aspergillus fumigatus' --target-taxid 330879 \
  --target-genome-id Afum --target-fasta "$WD/ref/Afum_Af293.fasta"

bash "$PIPELINE_DIR/rrmg_download_genomes.sh" \
  --manifest "$WD/input/rrmg/community_manifest.tsv" \
  --out-manifest "$WD/input/rrmg/community_manifest.downloaded.tsv"

# Inspect before continuing; successful background downloads are renormalized.
cat "$WD/input/rrmg/download_failed.tsv"
test -s "$WD/ref/Afum_Af293.fasta"

# Export only tool resolution settings for this direct utility invocation.
source "$PIPELINE_DIR/rrmg_environment.sh"
bash "$PIPELINE_DIR/rrmg_test_randomreadsmg_genomes.sh" \
  --manifest "$WD/input/rrmg/community_manifest.downloaded.tsv" \
  --out-manifest "$WD/input/rrmg/community_manifest.tested.tsv" \
  -w "$WD" --test-pairs 200 --readlen 150 --seed 999
```

The manifest columns are `genome_ID name taxid rel_abund fasta role` (TSV).
The target row has `role=target` and abundance 0; its spike level is set in the plan.
The Af293 example uses the [strain taxid 330879](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=330879).
A species-level taxid alone is not a promise of a particular strain. Keep
existing experimental IDs unchanged when resuming; audit reference identity separately.

The downloader uses [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_taxon/)
from PATH. It reuses non-empty FASTA files, tries reference genomes first, then
complete/chromosome/scaffold/contig packages. It chooses one assembly per row
from the package; fallback packages can contain many assemblies and be large.
`--threads` is currently a compatibility option, not download parallelism.
Diagnostics are saved in `input/rrmg/download_logs/`; failed target downloads
stop the command. Inspect failed background rows and the genome-test summary
before accepting the renormalized community. Preserve downloaded FASTA checksums
and NCBI assembly reports (`--keep-tmp`) to freeze references for a real study.
The target gene FASTA is **not** extracted by this downloader: provide your
validated `ref/CYP51A_wt.fasta` before Step3/4.

## 4. Create plans

These are **small pilot settings**, not the original 480-row design.
Use a separate work directory for the pilot; never overwrite an existing plan.

```bash
python3 "$PIPELINE_DIR/rrmg_make_plan.py" \
  --manifest "$WD/input/rrmg/community_manifest.tested.tsv" \
  --out "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  --sample-prefix Pilot_Afum \
  --target-name 'Aspergillus fumigatus' --target-taxid 330879 \
  --target-genome-id Afum --target-fasta "$WD/ref/Afum_Af293.fasta" \
  --target-abundances '0.01' --total-pairs '1000000' --seeds '1-2' --readlen 150

python3 "$PIPELINE_DIR/rrmg_make_background_plan.py" \
  --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  -w "$WD" --outdir "$WD/temp/rrmg_background" --background-seed 1000
```

This design yields one shared background condition and two permutation rows.
A low-depth pilot tests execution, not expected CYP51A recovery.
Record the actual FASTQ pair counts, tool versions and design parameters before
scaling up; the runner does not certify the simulator's read-count semantics.
Some RandomReadsMG releases interpret `reads=N` per contig rather than per
genome. For multi-contig targets compare Step5 `TargetPairsFastq` with
`TargetPairsDesign`; do not assume designed and realized abundances are equal.
This release does not silently change the simulation model for existing studies.

## 5. Run directly (no PBS)

Activate the environment once, then run stages in order. No extra for-loop is
needed for a contiguous range.

```bash
conda activate fungamr_meta
cd "$PIPELINE_DIR"
bash ./rrmg_submit.sh step1
bash ./rrmg_submit.sh step2 1-2
bash ./rrmg_submit.sh step3 1-2
bash ./rrmg_submit.sh step4 1-2
bash ./rrmg_submit.sh step5 1-2
# step45 1-2 is equivalent to step4 followed by step5 for each row.
# Optional, after QC and assembly: bash ./rrmg_submit.sh step3b 1-2
```

After `chmod +x rrmg_submit.sh`, `./rrmg_submit.sh step2 1-100` is identical.
`rrmg run step2 1-100 --config /path/config.sh` and `rrmg_run_local.sh` remain supported.
Selection is the original `run_index`, or **bg_index for step1**; it is not
Seed, target-species ID or a text line number including the header. A single
integer runs exactly one row. Omitted range runs all IDs in that stage's table.
Missing/duplicate IDs fail before launch; sparse tables retain original IDs.

```bash
bash ./rrmg_submit.sh step2 409-480 --dry-run
bash ./rrmg_submit.sh step2 409-480 --jobs 1
bash ./rrmg_submit.sh step2 225-230 --jobs 1
# Explicit replacement, only after inspecting the selected rows:
# bash ./rrmg_submit.sh step2 409-480 --force
```

Do not run overlapping invocations on the same stage/row or regenerate a
background while another job consumes it. Wait for upstream completion.
Dry-run validates the selection and dispatcher, not all FASTQ/reference inputs.
Per-invocation `logs/local/.../summary.tsv` records OK/FAILED and per-row logs.
Any row failure causes nonzero final exit; other selected rows are still attempted.
On an unscheduled remote machine use tmux/screen for long runs. On an HPC
cluster run analysis on allocated compute nodes, not the login node.

## 6. Ordinary PBS or PBS array

Copy and adapt [examples/pbs.pbs](examples/pbs.pbs) or
[examples/pbs_array.pbs](examples/pbs_array.pbs); these are not site-ready jobs.
They activate `fungamr_meta`, optionally select Java and call the same runner.

```bash
# Ordinary job: its body runs step2 1-100 sequentially.
qsub -v RRMG_STAGE=step2,RRMG_ROWS=1-100 examples/pbs.pbs
# Array-capable PBS Pro: each array task runs ONE original row.
qsub -J 1-100 -v RRMG_STAGE=step2 examples/pbs_array.pbs
```

In an array the command is `rrmg_submit.sh step2 "$PBS_ARRAY_INDEX"`,
**never** `step2 "1-$PBS_ARRAY_INDEX"`. Torque may use `-t` and `PBS_ARRAYID`;
the example accepts either variable. Edit queue/account/resource syntax for
your site; array limits/concurrency controls vary. Match ranges to bg_index
for step1, run_index for other stages. Root-level legacy `*_array.pbs` files
remain optional compatibility examples and are not used by the core.
Do not submit all stages together without validated dependencies.

## 7. Resume, collect and preserve results

Defaults reuse outputs according to the existing stage checks (often non-empty
files, not checksums). Step3 writes STEP3_SUCCESS after mapping/summary completion;
Step3B uses VAMB_SUCCESS. Step4/5 reuse per-run summaries. Keep target-only
reads (`STEP2_KEEP_TARGET_READS=1`) until Step5. Nothing in the new runner
automatically deletes reads. `--force` replaces selected outputs and can
remove assembly/binning directories; use only for inspected rows.

```bash
source "$RRMG_CONFIG"
python3 "$PIPELINE_DIR/rrmg_collect_map_results.py" \
  --summary "$WD/result/rrmg_split_pipeline_summary.tsv" \
  --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  --root "$WD/temp/rrmg_permutation" \
  --out "$WD/result/rrmg_target_gene_map_merged.tsv" --wd "$WD"
bash "$PIPELINE_DIR/rrmg_collect_supplementary_downstream.sh" \
  "$WD" "$WD/temp/rrmg_permutation"
# Only if VAMB was run:
# python3 "$PIPELINE_DIR/rrmg_collect_vamb_results.py" \
#   --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
#   --root "$WD/temp/rrmg_permutation" \
#   --out "$WD/result/rrmg_vamb_per_run.tsv" --wd "$WD"
```

Partial cohorts can be collected. Compare IDs and missing fields against the plan
before plotting; absent output is not biological zero. Back up depth/window files
as well as merged summaries if figures require them. No VAMB output is required
for a reads/contigs-only report.

## Testing and reproducibility

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
bash tests/smoke_test_local_runner.sh
python -m compileall -q .
find . -type f \( -name '*.sh' -o -name '*.pbs' \) -print0 | xargs -0 -n1 bash -n
python -m pip wheel --no-deps . -w dist
python tests/check_installed_wheel.py dist/*.whl
```

Tests use temporary synthetic fixtures/mocks, not real scheduler submissions.
Keep code, environment locks, small design templates and reference checksums
in version control; keep FASTA/FASTQ/BAM/results in managed data storage.
Pin a verified commit and do not update code/config during active runs.
Choose an appropriate license before distributing under an open-source license;
this update does not select one for the author.
