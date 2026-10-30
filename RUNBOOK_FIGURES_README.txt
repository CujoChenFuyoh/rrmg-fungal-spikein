RRMG figure inputs and rendering
===============================

Summary figures: rrmg_plot_preview.R
  Input: project plan, collected result TSVs, per-run gene depths.
  Output: WD/figure_output/local_preview/ (PNG/PDF, source data, sample audit).

Read/contig alignment figures: rrmg_plot_alignment.R
  Input: COMPLETE.json and plot_manifest.with_target_only.tsv from
         WD/figure_data/CYP51A_alignment/.
  Output: WD/figure_output/alignment_sites/ (depth/contig panels, read blocks,
          multi-page figures and optional reference-site coverage tables).

Three alignment sources are required for complete read/contig figures:
  contigs_to_GENE.paf                     contig -> target gene
  reads_to_GENE.paf.gz + depth.tsv         mixed-QC reads -> target gene
  target_only_reads_to_GENE.paf.gz + depth target-only reads -> target gene
Step5's whole-genome mapping does not replace target-only gene alignment.
Step4 now creates both mixed-read and target-only-read gene PAF/depth outputs.
Retain target-only FASTQs until Step4 succeeds; collection itself reads no FASTQ.
The full workflow keeps PAF output enabled; --no-paf is not its default.

Collect every plotting input (Linux server or Ubuntu; no mapping and no R):
  bash "$CODE/rrmg_collect_figure_data.sh" \
    --wd "$WD" --gene "$WD/ref/CYP51A_wt.fasta"
  Optional annotations: add --sites FILE.tsv --catalogue FILE.tsv together.
  For base-level evidence completeness, add --require-base-evidence.
  Step3/4 now add cs PAF and raw/Q20 base tables automatically; legacy projects
  can run rrmg_backfill_cyp51a_evidence.sh first. See docs/BASE_EVIDENCE.md.
  Existing R plots still show coverage, not validated SNV/FungAMR calls.
  This complete bundle also requires Step5 summary/window outputs for every run.
  The gene reference must match sequence orientation and reference coordinates.
  Invalid annotations are recorded in site_annotation_audit.tsv.
  Download the complete WD/figure_data directory to the local Windows WD.

Render with native R (use that R installation's packages):
  Rscript "$CODE/rrmg_plot_preview.R" "$WD" "$WD/figure_output/local_preview"
  Rscript "$CODE/rrmg_plot_alignment.R" "$WD/figure_data/CYP51A_alignment" \
    "$WD/figure_output/alignment_sites"
  Alignment dependencies: ggplot2, dplyr, tidyr, purrr, scales, patchwork.

Render both figure families in local Windows PowerShell:
  & "$CODE\rrmg_plot_all_windows.ps1" -WorkDir $WD -CodeDir $CODE
Multiple installations require -Rscript 'C:\Program Files\R\R-x.y.z\bin\Rscript.exe'.
Windows packages stay in Windows R. No package is installed by this launcher.

Outputs use fixed paths. Each helper stages all work, then replaces only its own
inventoried directory after success. Changed/added user files are protected.
Mapping/rendering failure keeps the previous successful output; require the
current invocation's success message, not an old COMPLETE.json.
Do not run concurrent writers into one managed directory.

Figure collection reuses audited Step3/Step4 outputs and never remaps reads or
contigs. Step4 resumes the mixed and target-only read sets independently, so an
older complete mixed-read mapping can be retained while the target-only mapping
is added.

Summary denominators show planned, evaluated and missing runs separately.
Missing data are not zero. Contig-display breadth and Step3 summary breadth use
the same existing Step3 PAF and MAPQ>=20 filter. Site labels show reference coverage,
not variant/phenotype calls. Read-block display samples up to350 alignments per
source; depth/site calculations use the complete depth data.
