# Workflow and output reference

## Stages

| Stage | Script | Input | Principal output |
| --- | --- | --- | --- |
| Manifest | `rrmg_prepare_manifest.py` | Taxon abundance table | Six-column manifest |
| Download | `rrmg_download_genomes.sh` | Manifest taxids | FASTA, updated manifest, failure list |
| Compatibility check | `rrmg_test_randomreadsmg_genomes.sh` | FASTA, manifest, BBMap | Tested manifest and logs |
| Plan | `rrmg_make_plan.py` | Tested manifest, target, abundance/depth/seeds | Stable `run_index` / `run_id` table |
| Background design | `rrmg_make_background_plan.py` | Plan, manifests, FASTA | Shared background plan/index/designs |
| Step1 | `rrmg_step1_make_background_reads.sh` | Background designs | Paired background FASTQ and pair-count audit |
| Step2 | `rrmg_step2_make_target_and_merge.sh` | Plan, background index/FASTQ, target genome | Exact mixed pairs, retained target pairs, generation audit |
| Step3 | `rrmg_step3_downstream_analysis.sh` | Mixed reads, gene FASTA | fastp QC, metaSPAdes contigs, contig-to-gene summary |
| Step4 | `rrmg_step4_reads_to_gene_depth.sh` | Mixed-QC reads, retained target-only reads, gene FASTA | Independent BAM, per-base depth, summary and PAF.gz for both read sets |
| Step5 (optional) | `rrmg_step5_afum_target_reads_genome_cov.sh` | Target-only reads, genome | Genome coverage and optional windows |
| Step45 | `rrmg_step4_5_supplementary_downstream.sh` | Step4/5 inputs | Both analyses; individual steps skip the other branch |
| Step3b (optional) | `rrmg_step3b_vamb_binning.sh` | QC reads, contigs, target reference | VAMB bins and recovery metrics |
| Collect | `rrmg_collect_*.py` | Existing outputs and plan | Merged TSVs |
| Summary figures | `rrmg_plot_preview.R` | Plan, summaries, per-base depths | Summary figures, source tables, missingness audit |
| Figure-data collection | `rrmg_collect_figure_data.sh` | Completed Step3-5 outputs, plan, gene | Portable summaries, three PAF sources, depths, bins, manifests and `COMPLETE.json` |
| Alignment figures | `rrmg_plot_alignment.R` | Complete alignment inputs | Contig/read blocks, depths and optional reference-site coverage |

Backgrounds are shared across seeds for one design condition. Step2 varies
target simulation/mixing seeds and records exact source/target/total counts.
FASTQ pairs, fastp consumed-read counts and mix-generation markers are checked
before downstream outputs are reused.

## Metrics and reference definitions

- Step3 contig-to-gene mapping uses minimap2 `asm5`; summary alignments require
  MAPQ>=20. `HasHit` follows the surviving alignments; it is not a predeclared
  clinical sensitivity/specificity threshold. No extra minimum gene breadth or
  identity criterion is imposed by that field.
- Step4 maps mixed-QC and target-only reads independently with
  `sr --secondary=no`, reports samtools-depth coverage and produces both PAFs
  for read-block display.
- Figure-data collection reuses Step3's contig PAF and applies the same
  MAPQ>=20 filter used by the Step3 breadth summary before local display.
- Target-only-to-gene alignment controls describe the simulated target source.
  Step5 aligns that source to the whole genome and answers a different question.
- The empty VCF / `VarN=0` in Step3 is a legacy placeholder for variant analysis
  that was not performed. Reference-site labels likewise show annotation/coverage.

## Mapping work and reuse

With complete gene mapping, each run invokes minimap2 five times: Step3=1,
Step4 mixed reads=2 and Step4 target-only reads=2. Step5 adds one whole-genome
call; optional species mapping/VAMB add other work. Figure-data collection adds
no minimap2 or samtools call.

The three PAF sources are required for complete read/contig block figures.
This workflow retains them by default. The low-level Step4 `--no-paf` option
only serves a coverage-only use case and is not the recommended figure workflow.
Step4 resumes the two read sets independently using the mix-generation audit.
Figure collection validates and copies the completed outputs; rerendering the
portable bundle does not require FASTQ, BAM, contigs or an aligner.

## Resource and completion rules

`rrmg run` is scheduler-free. The installed CLI and compatibility Bash launcher
share the same runner and dispatcher. `--jobs` is local concurrency, not an array.
The thread/memory options are per row and do not request scheduler resources.

Step3 includes QC. `--skip-qc` uses raw merged reads; it does not mean reuse QC.
Audited matching QC is reused by the normal Step3 path. Read-only operation from
scratch still requires a dedicated QC stage; `--skip-assemble` currently requires
existing contigs. `--force` deliberately replaces selected stage outputs.

Keep references and parameters fixed in an existing analysis. Current generation
checks are useful but do not fingerprint every changed reference or parameter.
For changed inputs, use an isolated project or explicitly verified rebuild;
do not treat an old completion marker as a reference-change audit.

Collection can cover a partial cohort. Check expected IDs and missing fields
against the plan; report evaluated denominators and keep missing distinct from zero.
