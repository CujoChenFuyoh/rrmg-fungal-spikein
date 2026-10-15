RRMG local runbook + Windows R figures
=====================================

Keep these files together:
  runbook.sh                         Sequential Ubuntu/WSL commands
  rrmg_windows_r.sh                  WSL -> Windows R, with Windows package library
  rrmg_plot_preview.R                Cohort/design summaries + gene depth previews
  rrmg_plot_common.R                 Safe fixed-directory output publication
  rrmg_prepare_alignment_inputs.py   Existing contigs/reads -> portable alignment inputs
  rrmg_plot_alignment.R              Existing ggplot/patchwork layout + reference-site overlays

Use runbook.sh interactively; do not run the whole file against an existing study.
Set paths in section 0. CODE is the repository or rrmg-exported core directory.
All helpers ship with the repository and pip package; no RUNBOOK_DIR is needed.
These helpers extend plotting only; they do not replace the core pipeline.

Step3 defaults to STEP3_ONLY_ASSEMBLER=1: fastp QC, then metaSPAdes without
BayesHammer. Set 0 in the active rrmg_config.sh to enable native correction.
The assembly mode is recorded; switching mode or reusing old unrecorded contigs
requires an explicit Step3 --force rebuild. Regenerate sections 11-15 afterward.
Do not treat file presence or a newly rendered figure as proof of a new Step3 run.

After Step1-5:
  11: collect results
  12: require complete IDs/files, exact target pairs, read-mix audits and fastp input counts
  13: Windows R summary figures
  14: prepare contig + mixed QC + target-only alignments, validate site annotations
  15: Windows R annotation/alignment-block figures and per-site summary tables

Windows R and packages:
  Rscript.exe: C:/Program Files/R/R-4.6.0/bin/Rscript.exe
  Package library: C:/Users/q6321/AppData/Local/R/win-library/4.6
Override WIN_RSCRIPT (WSL path) and WIN_R_LIB (Windows path) on another workstation.
The wrapper converts filesystem arguments with wslpath, including Ubuntu UNC paths.
No global R configuration or package installation is changed.

Input annotations:
  CYP51A_FungAMR_sites.tsv and Afum_fungamr.catalogue.tsv in the configured DB directory.
  NCBI chromosome interval in the gene FASTA header must match its orientation and length.
  Exact codon coordinates and reference bases are checked. Excluded entries are reported.
  This is reference-site annotation/coverage, not mutation calling or an AMR phenotype.

Outputs:
  $WD/figure_output/local_preview/: PNG, PDF, source TSVs, sample audit, checksums.
  $WD/figure_data/CYP51A_alignment/: portable FASTA/PAFs/depths/manifest/site audit.
  $WD/figure_output/alignment_sites/: PNG/PDF, multi-page PDFs and per-site tables.
No dates are appended. A successful rerun replaces only the inventoried outputs
of the same helper. New rendering/mapping completes in staging before publication.
Old extra pages/rows are removed. Unmanaged directories or added/modified files
are protected: preserve/move that directory first rather than deleting its guard.
Failures keep the previous successful outputs; require the current [DONE] message.
Do not run simultaneous writers for one output directory. After an interrupted
process, check that no writer remains before removing its empty .lock directory.
Section 14 rebuilds the selected rows, not an incremental append. For replotting
unchanged inputs, skip section 14 and run section 15 only.
Retain Step2 target-only FASTQs and Step3 QC FASTQs until section 14 finishes.
Original contigs, FASTQs, BAMs, annotation tables and old scripts are not deleted.

Section 14 read-block inputs (row_00001/, row_00002/, ...):
  contigs_to_GENE.paf
  reads_to_GENE.paf.gz                  mixed QC reads mapped to the gene
  target_only_reads_to_GENE.paf.gz      target-only reads mapped to the gene
  reads_to_GENE.depth.tsv and target_only_reads_to_GENE.depth.tsv
COMPLETE.json and plot_manifest.with_target_only.tsv appear after all rows succeed.
Step5 target-genome alignments are not a substitute for the target-only gene PAF.
An IsADirectoryError for '.' before [MAP] usually means an empty shell variable;
set CATALOGUE_SOURCE to the actual Afum_fungamr.catalogue.tsv, not its directory.

Summary denominators are evaluable runs, with planned and missing counts shown.
Zero is not substituted for missing data. Fractions from a small pilot are descriptive.
Read-block displays are limited to 350 alignments per source; per-site depth uses all bases.
Alignment-span coverage does not prove every base or codon is correct.

Validation on 2026-10-14:
  Bash syntax checks passed.
  Eight Python reference-coordinate/input/output tests passed (29 repository tests total).
  Wheel installation/export and source-distribution rebuild include all plotting assets.
  Two consecutive fixed-directory mapping and Windows-R plotting runs passed.
  Empty catalogue variables, user-added/modified outputs and failed reruns were tested.
  R export/ID/denominator/missingness tests passed on Windows R and Ubuntu R.
  WSL -> Windows R export passed from /mnt/c and the Ubuntu filesystem via UNC.
  Two existing random-DNA simulation replicates completed new alignment preparation.
  Artificial TEST_SITE_A/B annotations exercised labelled figures and per-site summaries.
  Existing real CYP51A annotations were validated separately: 50 included, 26 excluded.
  PNGs were visually checked and label/axis clipping corrected.
This is software validation, not a rerun or performance assessment of the biological study.
