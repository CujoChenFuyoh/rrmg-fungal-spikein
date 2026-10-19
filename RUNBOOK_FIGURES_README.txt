RRMG figure inputs and rendering
===============================

Summary figures: rrmg_plot_preview.R
  Input: project plan, collected result TSVs, per-run gene depths.
  Output: WD/figure_output/local_preview/ (PNG/PDF, source data, sample audit).

Read/contig alignment figures: rrmg_plot_alignment.R
  Input: COMPLETE.json and plot_manifest.with_target_only.tsv from
         rrmg_prepare_alignment_inputs.py.
  Output: WD/figure_output/alignment_sites/ (depth/contig panels, read blocks,
          multi-page figures and optional reference-site coverage tables).

Three alignment sources are required for complete read/contig figures:
  contigs_to_GENE.paf                     contig -> target gene
  reads_to_GENE.paf.gz + depth.tsv         mixed-QC reads -> target gene
  target_only_reads_to_GENE.paf.gz + depth target-only reads -> target gene
Step5's whole-genome mapping does not replace target-only gene alignment.
Retain PAFs, QC reads and target-only reads until figure preparation succeeds.
The full workflow keeps PAF output enabled; --no-paf is not its default.

Prepare alignment inputs (active Linux analysis environment):
  python3 "$CODE/rrmg_prepare_alignment_inputs.py" \
    --wd "$WD" --gene "$WD/ref/CYP51A_wt.fasta" \
    --out "$WD/figure_data/CYP51A_alignment" --rows 1-2 --threads 4
  Optional annotations: add --sites FILE.tsv --catalogue FILE.tsv together.
  The gene reference must match sequence orientation and reference coordinates.
  Invalid annotations are recorded in site_annotation_audit.tsv.

Render with native R (use that R installation's packages):
  Rscript "$CODE/rrmg_plot_preview.R" "$WD" "$WD/figure_output/local_preview"
  Rscript "$CODE/rrmg_plot_alignment.R" "$WD/figure_data/CYP51A_alignment" \
    "$WD/figure_output/alignment_sites"
  Alignment dependencies: ggplot2, dplyr, tidyr, purrr, scales, patchwork.

Or render with Windows R from WSL:
  export WIN_RSCRIPT='/mnt/c/Program Files/R/R-<version>/bin/Rscript.exe'
  # Optional: WIN_R_LIB is a Windows library path; otherwise use R's own defaults.
  bash "$CODE/rrmg_windows_r.sh" --check
  bash "$CODE/rrmg_windows_r.sh" "$CODE/rrmg_plot_preview.R" \
    "$WD" "$WD/figure_output/local_preview"
  bash "$CODE/rrmg_windows_r.sh" "$CODE/rrmg_plot_alignment.R" \
    "$WD/figure_data/CYP51A_alignment" "$WD/figure_output/alignment_sites"
The launcher detects Rscript.exe on PATH or a single standard Windows installation.
Multiple installations require WIN_RSCRIPT. Windows packages stay in Windows R;
Linux R does not load Windows compiled packages. No package is installed by this launcher.

Outputs use fixed paths. Each helper stages all work, then replaces only its own
inventoried directory after success. Changed/added user files are protected.
Mapping/rendering failure keeps the previous successful output; require the
current invocation's success message, not an old COMPLETE.json.
Do not run concurrent writers into one managed directory.

Figure preparation currently remaps its selected cohort; it is not incremental.
Once preparation succeeds, rerendering the same inputs does not require mapping.
Reference/generation/parameter-aware reuse is tracked in docs/REVIEW_CHECKLIST.zh-CN.md.

Summary denominators show planned, evaluated and missing runs separately.
Missing data are not zero. Contig-display breadth and Step3 summary breadth use
different filters; see docs/WORKFLOW.md. Site labels show reference coverage,
not variant/phenotype calls. Read-block display samples up to350 alignments per
source; depth/site calculations use the complete depth data.
