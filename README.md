# RRMG fungal spike-in workflow

Version `0.3.3`; release notes are in [CHANGELOG.md](CHANGELOG.md).

RRMG evaluates target-gene recovery across metagenomic background, target
read fraction, sequencing depth and random-seed combinations. It generates
audited paired reads, assembles contigs, measures gene coverage and exports
summary and alignment-block figures.

## Start here

| Task | Guide |
| --- | --- |
| First run / smallest set of commands | [Quick start (中文)](docs/QUICKSTART.zh-CN.md) |
| Environment, manifest, genome download and plan | [Installation and inputs (中文)](docs/INSTALL.zh-CN.md) |
| Stage inputs, outputs and metric definitions | [Workflow reference](docs/WORKFLOW.md) |
| Interactive pilot from inputs through figures | [runbook.sh](runbook.sh) |
| Figure inputs, native R and Windows R in WSL | [Figure guide](RUNBOOK_FIGURES_README.txt) |
| Ordinary PBS / PBS array | [Scheduler examples](pbs_example/README.md) |
| Failed runs, old outputs and provenance | [Troubleshooting (中文)](docs/TROUBLESHOOTING.zh-CN.md) |
| Verified redundancy and changes awaiting confirmation | [Review checklist (中文)](docs/REVIEW_CHECKLIST.zh-CN.md) |

## One analysis entry point

Activate the analysis environment and install the checked-out source:

```bash
conda activate fungamr_meta
python -m pip install .
rrmg init --workdir "$HOME/SpikeIN" --pipeline-dir "$PWD"
export RRMG_CONFIG="$HOME/SpikeIN/rrmg_config.sh"
rrmg doctor
```

`pip` installs workflow code. Bioinformatics programs and experimental inputs
are prepared separately, as described in the installation guide. `init` creates
one configuration and preserves existing files. Run it once per new project.

After preparing a plan and references:

```bash
rrmg run step1 -t 4 -m 4 --jobs 1
rrmg run step2 1-2 -t 1 -m 4 --jobs 1
rrmg run step3 1-2 -t 4 -m 16 --jobs 1
rrmg run step4 1-2 -t 4 --jobs 1
# Required by the complete figure-data bundle; optional for gene-only analysis:
rrmg run step5 1-2 -t 4 --jobs 1
```

Step3 includes fastp QC, metaSPAdes and contig-to-gene mapping. Step4 maps both
mixed metagenome QC reads and retained target-only reads to the target gene,
retaining PAF and per-base depth for alignment-block plots. Step5 remains the
separate target-only-read to whole-genome coverage analysis.

`rrmg run` executes inside the current allocation and never submits jobs.
`--jobs` controls concurrent rows; `-t` and `-m` apply per row. Put the same
command in a site-specific PBS/Slurm shell to schedule it. On an HPC, run heavy
analysis on compute nodes. `bash rrmg_submit.sh ...` is the source-checkout
compatibility entry; `rrmg_run_local.sh` is its implementation, not a separate
analysis engine. New scheduler examples live only in `pbs_example/`.

## Choose the outputs you need

- **Assembly recovery:** Steps1–3, then collect contig mapping results.
- **Read and contig comparison:** add Step4 and retain PAF/depth outputs.
- **Alignment/site figures:** Step3 plus Step4 produce all three gene-alignment
  sources. Collect a portable `figure_data/` bundle, then render it locally.
- **Additional diagnostics:** Step5 whole-genome coverage, Kraken2 and VAMB
  (`step3b`) are optional. Step5 coverage windows are not metagenomic bins.

The unified all-figure collector includes the Step5 whole-genome panels and
therefore requires completed Step5 summaries/windows. Step5 remains optional if
only the gene-level Step3/4 outputs are needed.

For complete read/contig figures, retain target-only reads through Step4. After
the complete `figure_data/` bundle has been verified and backed up, plotting no
longer needs FASTQ, BAM or contig FASTA files.

## Data layout and interpretation

Code and analysis data are separate. A project uses one `rrmg_config.sh`,
`input/rrmg/`, `ref/`, `temp/rrmg_background/`, per-run directories under
`temp/rrmg_permutation/`, `result/`, `figure_data/`, `figure_output/` and `logs/`.
Existing projects may configure another background directory.

The supported abundance mode is `read_fraction`. Pair-count, fastp-input and
mix-generation checks protect comparisons between consistent outputs. Keep
reference sequences and method parameters fixed during a run; reference/parameter
fingerprint improvements remain listed in the review checklist.

Gene recovery and reference-site coverage are sequence-based measurements.
Step3's current empty VCF and `VarN=0` are placeholders for an unperformed
variant analysis. Partial cohorts can be collected; missing outputs remain
missing rather than biological zeros. A successful collector is not a cohort
completion check. Metric definitions and figure denominators are documented in
the workflow and figure guides.

## Development and reproducibility

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
bash tests/smoke_test_local_runner.sh
python -m pip wheel --no-deps . -w dist
python tests/check_installed_wheel.py dist/*.whl
```

Tests use small fixtures and mocked tools unless explicitly stated. Pin the
code commit, reference checksums and environment; retain validation records
with those inputs. Keep sequencing data and private site configuration outside
the repository. Do not upgrade code/configuration during an active analysis.

## License

This repository currently has no `LICENSE` file. A project license remains
to be selected by the author.
