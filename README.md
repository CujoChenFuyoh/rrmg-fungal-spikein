# RRMG fungal spike-in permutation workflow

Version: `0.2.0`

Repository: [CujoChenFuyoh/rrmg-fungal-spikein](https://github.com/CujoChenFuyoh/rrmg-fungal-spikein).

This repository simulates paired-end metagenomic reads, inserts a fungal target
at controlled relative abundances, assembles each permutation, measures target
gene recovery, and optionally performs one-sample VAMB binning. The reference
implementation uses *Aspergillus fumigatus* and CYP51A, but the plan-generation
scripts accept another target genome and target gene.

The maintained workflow is the split Step 1–5 implementation. It reuses one
background read set for every abundance × sequencing-depth condition and varies
only the target reads among permutation seeds.

> **Current limitation:** Step 3 measures target-gene presence, alignment
> identity, and breadth. It writes an empty VCF placeholder and reports
> `VarN=0`; it does **not** call or annotate fungAMR resistance mutations.
> Mutation calling must be added as a separately validated stage before these
> results can be interpreted as AMR genotypes.

## Workflow

| Stage | Purpose | Main program | Input | Principal output |
|---|---|---|---|---|
| Prepare | Convert a taxon table to a community manifest | `rrmg_prepare_manifest.py` | Taxon abundance table | `community_manifest.tsv` |
| Prepare | Download and test reference genomes | `rrmg_download_genomes.sh`, `rrmg_test_randomreadsmg_genomes.sh` | Manifest | Tested manifest and FASTA files |
| Plan | Expand abundance × depth × seed combinations | `rrmg_make_plan.py` | Tested manifest and target FASTA | `rrmg_permutation_plan.tsv` |
| Plan | Collapse seed replicates into reusable backgrounds | `rrmg_make_background_plan.py` | Permutation plan | `background_plan.tsv`, `background_index.tsv` |
| Step 1 | Simulate one background read set per condition | `rrmg_step1_make_background_reads.sh` | Background plan and genomes | Background FASTQ pair |
| Step 2 | Simulate target reads and concatenate with background | `rrmg_step2_make_target_and_merge.sh` | Plan, background FASTQ, target FASTA | Mixed FASTQ pair and merge design |
| Step 3 | QC, optional Kraken2, metaSPAdes, contig-to-gene mapping | `rrmg_step3_downstream_analysis.sh` | Mixed FASTQ and target gene | QC reads, contigs, PAF, per-run summary |
| Step 3B | Bin each permutation independently with VAMB | `rrmg_step3b_vamb_binning.sh` | QC reads, contigs, target genome | VAMB bins and recovery metrics |
| Step 4/5 | Optional read-to-gene depth and target-genome coverage | `rrmg_step4_5_supplementary_downstream.sh` | Mixed/target-only reads | Depth and coverage summaries |
| Collect | Merge all per-run results | `rrmg_collect_*.py`, `rrmg_collect_supplementary_downstream.sh` | Plan and per-run outputs | Combined TSV files |

Each permutation keeps the existing directory layout:

```text
WORK_DIR/
├── input/rrmg/
├── ref/
├── temp/
│   ├── rrmg_background/
│   └── rrmg_permutation/<run_id>/
│       ├── design/
│       ├── reads/
│       ├── qc/
│       ├── kraken/
│       ├── asm/contigs.fasta
│       ├── map/
│       └── binning/
├── result/
└── logs/pbs/
```

## Installation

The workflow targets Linux, GNU Bash, and Python >=3.10. PBS/PBS Pro is optional;
the same original scripts can run through the local runner without job arrays.
The `rrmg` command is a thin installation/dispatch layer, not a new analysis engine.

**Chinese migration/installation guide:** [docs/INSTALL.zh-CN.md](docs/INSTALL.zh-CN.md).

Activate your existing analysis environment, clone the repository, and install:

```bash
conda activate fungamr_meta
git clone https://github.com/CujoChenFuyoh/rrmg-fungal-spikein.git rrmg-fungal-spikein-github
cd rrmg-fungal-spikein-github
python -m pip install .
rrmg --version
```

The explicit checkout name `rrmg-fungal-spikein-github` preserves the directory
used in the existing UBDA configuration; it need not match the GitHub repository
name. If this checkout already exists, use it instead of cloning over it.

Alternatively, install directly from the repository's `main` branch:

```bash
python -m pip install "git+https://github.com/CujoChenFuyoh/rrmg-fungal-spikein.git@main"
```

For reproducible runs, pin a tested commit that has been pushed to this repository.
Run the following from that commit's local checkout:

```bash
COMMIT_SHA="$(git rev-parse HEAD)"
python -m pip install "git+https://github.com/CujoChenFuyoh/rrmg-fungal-spikein.git@${COMMIT_SHA}"
```

Git and repository access are required for Git URL installs. For a private
repository, authenticate Git on the target machine first; a browser login does
not automatically authenticate command-line Git. Do not put tokens in these URLs.
Local, uncommitted edits are only installed with `pip install .`, not by a Git URL.

This project has not automatically been published to PyPI.
pip installs the core code, **not** external bioinformatics programs or
input/reference/result data. If you need a new analysis environment, create it
from the cloned repository first:

```bash
conda env create -f environment.yml
conda activate fungamr_meta
python -m pip install .
```

`environment.yml` is intentionally human-readable. For a frozen Linux build,
export the exact tested package set and commit it with the release:

```bash
conda list -n fungamr_meta --explicit > conda-linux-64.lock
# Recreate it later with:
conda create -n fungamr_meta --file conda-linux-64.lock
```

## Configuration

Existing project/configuration? Keep it; no initialization or plan regeneration is needed:

```bash
export RRMG_CONFIG=/path/to/rrmg_config.sh
rrmg run step2 451-455 --jobs 1 --dry-run
# Run analysis only on allocated compute resources:
rrmg run step2 451-455 --jobs 1
```

For a new site config bound to a cloned repository:

```bash
rrmg init --workdir /path/to/SpikeIN --pipeline-dir "$PWD" --env fungamr_meta
rrmg doctor --workdir /path/to/SpikeIN
```

Without `--pipeline-dir`, `init` exports the original core scripts to a versioned
`.rrmg/core-0.2.0/` directory inside the project and pins the new config to it.
No input data or existing config is replaced; unresolved software paths are
marked `/EDIT_ME/`. Review paths and resources before analysis.
`rrmg export /path/to/new_core` copies core scripts and a SHA-256 manifest to a new
directory. `rrmg path` prints the installed script location. Neither command
copies sequencing data or installs external tools.

The original manual configuration and Bash entry points remain supported:

Machine-specific paths are not stored in PBS files. Create one local config:

```bash
cp rrmg_config.example.sh rrmg_config.sh
nano rrmg_config.sh
```

Set at least:

- `PIPELINE_DIR`: this cloned repository;
- `WD`: working directory containing `input/`, `ref/`, `temp/`, and `result/`;
- `CONDA_SH` and `RRMG_CONDA_ENV`;
- `RANDOMREADSMG_SH`, `JAVA_HOME`, and `JAVA_BIN`;
- reference and Kraken2 paths as required.

The private `rrmg_config.sh` is ignored by Git. Verify the installation before
submitting expensive jobs:

```bash
bash rrmg_preflight.sh
```

The `#PBS -l` resource declarations remain in each PBS wrapper because PBS
reads directives before Bash can load the config. If thread or memory values
are changed in `rrmg_config.sh`, submit matching PBS resources with `qsub -l`
or update the wrapper deliberately.

## Required inputs

The canonical community manifest has six tab-separated columns:

```text
genome_ID  name  taxid  rel_abund  fasta  role
```

`role` is `background` or `target`. Keep the target row in the manifest with
`rel_abund=0`; the actual target abundance is specified when creating the
permutation plan.

At minimum, provide:

- a tested community manifest and all referenced genome FASTA files;
- the target genome FASTA;
- the target-gene FASTA;
- a Kraken2 database only when `RUN_KRAKEN=1`.

For a new experimental design, generate a plan from recorded parameters so every
run has a unique `run_index` and `run_id`. For migration/resumption, retain the
existing plan and identities; audit old absolute paths instead of regenerating
the experiment merely because its code was installed on a new server.

## Create plans

Run commands from the code repository after loading the local configuration:

```bash
source rrmg_config.sh
source "$CONDA_SH"
conda activate "$RRMG_CONDA_ENV"
mkdir -p "$WD/input/rrmg" "$WD/temp" "$WD/result" "$WD/logs/pbs"

python3 rrmg_make_plan.py \
  --manifest "$WD/input/rrmg/community_manifest.tested.tsv" \
  --out "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  --sample-prefix AirMG_Afum \
  --target-name "Aspergillus fumigatus" \
  --target-taxid 746128 \
  --target-genome-id Afum \
  --target-fasta "$WD/ref/Afum_Af293.fasta" \
  --target-abundances "0.005,0.01,0.02" \
  --total-pairs "10000000,25000000,50000000" \
  --seeds "1-100" \
  --readlen 150

python3 rrmg_make_background_plan.py \
  --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  -w "$WD" \
  --outdir "$WD/temp/rrmg_background" \
  --background-seed 1000
```

The numerical values above are examples. Record the values used for each
analysis in a release, methods document, or checked-in small design file.

## Run with PBS

`rrmg_submit.sh` calculates array bounds from the current plan, so stale
hard-coded `#PBS -J` ranges are not used:

```bash
./rrmg_submit.sh step1
# Wait for Step 1 to finish successfully.
./rrmg_submit.sh step2
# Wait for Step 2.
./rrmg_submit.sh step3
# Wait for Step 3 assemblies.
./rrmg_submit.sh step3b

# Optional supplementary analyses:
./rrmg_submit.sh step45
```

Submit a selected range when recovering failed rows:

```bash
./rrmg_submit.sh step2 215-480
```

Extra PBS options are passed through after the optional range:

```bash
./rrmg_submit.sh step3 1-300 -W depend=afterok:12345
```

Do not submit a downstream stage until its required upstream outputs are
complete. Job IDs and plan rows should be recorded with the analysis metadata.

## Run without PBS

`rrmg_run_local.sh` runs the same maintained stage wrappers directly and does
not call `qsub`. Run one stage at a time and wait for it to finish before
starting its downstream stage:

```bash
bash rrmg_run_local.sh step1
bash rrmg_run_local.sh step2
bash rrmg_run_local.sh step3
bash rrmg_run_local.sh step3b

# Optional supplementary analyses:
bash rrmg_run_local.sh step45
```

Select one plan row or a range in the same way as for PBS arrays:

```bash
bash rrmg_run_local.sh step2 215
bash rrmg_run_local.sh step2 215-480
```

The default is one row at a time. Limited parallel execution is available,
but every simultaneous row uses the full per-row CPU and memory settings:

```bash
bash rrmg_run_local.sh step2 215-480 --jobs 4
bash rrmg_run_local.sh step3 215-220 --jobs 1
```

Use `--dry-run` to validate the range without starting analysis. Use `--force`
only after inspecting the selected rows; it overrides the corresponding force
setting from `rrmg_config.sh` for that invocation:

```bash
bash rrmg_run_local.sh step2 215-220 --dry-run
bash rrmg_run_local.sh step2 215-220 --force
```

Each invocation writes a status table and per-row logs under
`$WD/logs/local/<run_id>/`. The command returns a non-zero exit status when any
row fails and prints the failed row numbers. On a server without a scheduler,
run long commands inside `tmux` or `screen` so they survive an SSH disconnect.

## Resume and replacement behavior

- Step 1 skips a background when both background FASTQ files are non-empty.
- Step 2 skips a run when both mixed FASTQ files are non-empty.
- Step 3 reuses existing QC, Kraken2, and assembly outputs and writes
  `STEP3_SUCCESS` only after mapping and summary update finish.
- Step 3B skips a run only when both `VAMB_SUCCESS` and its summary exist.
- Step 4/5 skip completed per-run summaries.
- Force flags are disabled in `rrmg_config.example.sh`. Enable a force setting
  only for inspected rows because replacement can remove an incomplete
  assembly or VAMB directory.

`--skip-assemble` means reuse a non-empty existing `asm/contigs.fasta`; it never
creates or truncates that file.

## Collect results

After all relevant arrays finish:

```bash
source rrmg_config.sh

python3 "$PIPELINE_DIR/rrmg_collect_map_results.py" \
  --summary "$WD/result/rrmg_split_pipeline_summary.tsv" \
  --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  --root "$WD/temp/rrmg_permutation" \
  --out "$WD/result/rrmg_target_gene_map_merged.tsv" \
  --wd "$WD"

python3 "$PIPELINE_DIR/rrmg_collect_vamb_results.py" \
  --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  --root "$WD/temp/rrmg_permutation" \
  --out "$WD/result/rrmg_vamb_per_run.tsv" \
  --wd "$WD"

bash "$PIPELINE_DIR/rrmg_collect_supplementary_downstream.sh" \
  "$WD" "$WD/temp/rrmg_permutation"
```

## Reproducibility and data safety

- Commit code, small manifests, plan-generation parameters, configuration
  examples, environment specifications, and checksums—not FASTQ/BAM/FASTA data.
- Keep reference databases and large outputs in managed storage and record
  their version, acquisition date, and SHA-256 checksum.
- Tag every tested release, for example `v0.1.0`; do not silently edit a tag.
- Save the exact Conda explicit lock from the validated Linux environment.
- Back up final assemblies and result tables outside temporary/ephemeral space.
- Run the included GitHub Actions syntax checks on every push and pull request.

Local syntax checks:

```bash
python3 -m compileall -q .
find . -type f \( -name '*.sh' -o -name '*.pbs' \) -print0 \
  | xargs -0 -n1 bash -n

# Lightweight dispatcher tests; no bioinformatics programs or real qsub run.
bash tests/smoke_test_local_runner.sh
```

## Repository scope

Historical code snapshots, the old monolithic runner, compatibility Step 4
wrappers, runtime results, references, and local configuration are retained on
the author's workstation but excluded from Git by `.gitignore`. The canonical
public interface is the split workflow documented here.

Before making the repository public, select an appropriate software license
(for example MIT, BSD-3-Clause, or GPL-3.0) and add author/contact and citation
information. No license is chosen automatically by this repository.
