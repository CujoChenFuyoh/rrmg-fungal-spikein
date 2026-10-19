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
| Step4 | `rrmg_step4_reads_to_gene_depth.sh` | QC reads, gene FASTA | BAM, per-base depth, coverage summary, read PAF.gz |
| Step5 (optional) | `rrmg_step5_afum_target_reads_genome_cov.sh` | Target-only reads, genome | Genome coverage and optional windows |
| Step45 | `rrmg_step4_5_supplementary_downstream.sh` | Step4/5 inputs | Both analyses; individual steps skip the other branch |
| Step3b (optional) | `rrmg_step3b_vamb_binning.sh` | QC reads, contigs, target reference | VAMB bins and recovery metrics |
| Collect | `rrmg_collect_*.py` | Existing outputs and plan | Merged TSVs |
| Summary figures | `rrmg_plot_preview.R` | Plan, summaries, per-base depths | Summary figures, source tables, missingness audit |
| Alignment inputs | `rrmg_prepare_alignment_inputs.py` | Contigs, QC and target-only reads, gene | Three PAF sources, depths, manifest, `COMPLETE.json` |
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
- Step4 maps mixed-QC reads with `sr --secondary=no`, reports the current
  samtools-depth coverage and produces a PAF for read-block display.
- Alignment-input preparation uses `asm5 --secondary=no -c` for contigs and
  does not apply Step3's MAPQ>=20 breadth filter. Its display breadth can therefore
  differ from the Step3 summary. Preserve these labels until a common definition
  is chosen; do not copy a contig PAF between stages as if definitions matched.
- Target-only-to-gene alignment controls describe the simulated target source.
  Step5 aligns that source to the whole genome and answers a different question.
- The empty VCF / `VarN=0` in Step3 is a legacy placeholder for variant analysis
  that was not performed. Reference-site labels likewise show annotation/coverage.

## Mapping work and reuse

With gene mapping, default Step4 PAF and full alignment preparation, each run
usually invokes minimap2 eight times: Step3=1, Step4=2, figure preparation=5.
Step5 adds one whole-genome call; optional species mapping/VAMB add other work.
The mixed-QC pair is scanned four times across Step4 and figure preparation.

The three PAF sources are required for complete read/contig block figures.
This workflow retains them by default. The low-level Step4 `--no-paf` option
only serves a coverage-only use case and is not the recommended figure workflow.
Currently figure preparation recomputes alignments; rerendering existing complete
figure inputs does not. Automatic reuse awaits matching fingerprints for reference,
read generation, parameters and tool versions; see the review checklist.

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
