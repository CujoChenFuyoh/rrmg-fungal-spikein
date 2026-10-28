#!/usr/bin/env Rscript
# Base-R regression checks for the local preview helper. No analysis tools needed.
# Usage: Rscript tests/test_rrmg_plot_preview.R /path/to/rrmg_plot_preview.R
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
source(args[1])

run_tests <- function() {
  root <- tempfile("rrmg-figure-test-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "input/rrmg"), recursive = TRUE)
  dir.create(file.path(root, "result"))
  dir.create(file.path(root, "figure_data/CYP51A_alignment/row_00017"), recursive = TRUE)
  plan <- data.frame(run_id = c("a", "b"), run_index = c("17", "83"), seed = c(1, 2),
                      target_abundance = c(0.01, 0.02), total_pairs = c(10000, 10000))
  write_tsv(plan, file.path(root, "input/rrmg/rrmg_permutation_plan.tsv"))
  contigs <- data.frame(RunID = c("b", "a"), TargetBreadthPct = c(60, 0), HasGeneHit = c("Y", "N"))
  contig_path <- file.path(root, "result/rrmg_target_gene_map_merged.tsv")
  write_tsv(contigs, contig_path)
  reads <- data.frame(RunID = "a", BreadthPct = 0, MeanDepth = 0, TargetGeneBp = 3, MappedReadAlignments = 0,
                       DepthTSV = "/old/nonexistent/server/path/depth.tsv")
  write_tsv(reads, file.path(root, "result/rrmg_reads_to_gene_depth_summary.tsv"))
  target <- data.frame(RunID = c("b", "a"), ObservedBreadthPct = c(90, 80), ObservedMeanDepth = c(2, 1),
                       TargetPairsDesign = c(200, 100), TargetPairsFastq = c(400, 100))
  write_tsv(target, file.path(root, "result/rrmg_afum_target_reads_genome_coverage_summary.tsv"))
  depth_path <- file.path(root, "figure_data/CYP51A_alignment/row_00017/reads_to_GENE.depth.tsv")
  writeLines(c("gene\t1\t0", "gene\t2\t0", "gene\t3\t0"), depth_path)
  first <- collect_preview(root)
  stopifnot(identical(first$data$RunID, c("a", "b")),
            identical(first$data$ContigGeneBreadthPct, c(0, 60)),
            first$data$ReadGeneMeanDepth[1] == 0, is.na(first$data$ReadGeneMeanDepth[2]),
            first$status$reads_summary[2] == "missing", first$status$gene_depth[2] == "missing",
            identical(first$status$target_pairs, c("match", "MISMATCH")),
            nrow(first$depth) == 3, first$status$depth_vs_summary[1] == "match")
  cat("[PASS] Sparse row IDs, reordered joins, real zero versus missing, relocated depth, pair mismatch.\n")
  rates <- design_summary(first$data)
  stopifnot(identical(rates$PlannedN, c(1L, 1L)),
            identical(rates$ReadValidN, c(1L, 0L)),
            rates$ReadHitPct[1] == 0, is.na(rates$ReadHitPct[2]),
            identical(rates$ContigHitPct, c(0, 100)),
            identical(rates$ReadMissingN, c(0L, 1L)))
  cat("[PASS] Design rates use evaluable denominators, with planned and missing counts retained.\n")

  bad <- rbind(contigs, contigs[1, ])
  write_tsv(bad, contig_path)
  stopifnot(inherits(try(collect_preview(root), silent = TRUE), "try-error"))
  write_tsv(contigs, contig_path)
  bad <- contigs
  bad$TargetBreadthPct[1] <- 101
  write_tsv(bad, contig_path)
  stopifnot(inherits(try(collect_preview(root), silent = TRUE), "try-error"))
  write_tsv(contigs, contig_path)
  cat("[PASS] Duplicate RunID and invalid percentage rejected.\n")

  writeLines(c("gene\t1\t0", "gene\t3\t0"), depth_path)
  gaps <- collect_preview(root)
  stopifnot(gaps$status$gene_depth[1] == "gaps", nrow(gaps$depth) == 2,
            gaps$status$depth_vs_summary[1] == "MISMATCH")
  writeLines(c("gene\t1\t0", "gene\t1\t0"), depth_path)
  invalid <- collect_preview(root)
  stopifnot(invalid$status$gene_depth[1] == "invalid", nrow(invalid$depth) == 0)
  writeLines(c("geneA\t1\t0", "geneA\t2\t0", "geneB\t1\t0"), depth_path)
  multi <- collect_preview(root)
  stopifnot(length(multi$panels) == 3L,
            identical(vapply(multi$panels[1:2], function(p) p$ref, character(1)), c("geneA", "geneB")))
  cat("[PASS] Depth gaps preserved, malformed data flagged, reference sequences separated.\n")

  out <- file.path(root, "figures")
  main(c(root, out))
  expected <- c("overview_01.png", "overview_01.pdf", "gene_depth_01.png", "gene_depth_01.pdf",
                 "summary_source_data.tsv", "gene_depth_source_data.tsv", "sample_status.tsv", "input_checksums.tsv",
                 "design_hit_rates.png", "contig_metrics_by_design_01.pdf", "read_metrics_by_design_01.pdf",
                 "target_genome_metrics_by_design_01.png", "target_expected_vs_observed.png", "detection_by_design.tsv")
  stopifnot(all(file.info(file.path(out, expected))$size > 0))
  saved <- read_tsv(file.path(out, "summary_source_data.tsv"))
  stopifnot(saved$ReadGeneMeanDepth[1] == "0", is.na(saved$ReadGeneMeanDepth[2]))
  main(c(root, out))
  stopifnot(all(file.info(file.path(out, expected))$size > 0))
  writeLines("user notes", file.path(out, "notes.txt"))
  stopifnot(inherits(try(main(c(root, out)), silent = TRUE), "try-error"))
  stopifnot(readLines(file.path(out, "notes.txt")) == "user notes")
  cat("[PASS] PNG/PDF/source exports, fixed-directory rerun, zero versus NA, added files protected.\n")

  managed <- file.path(root, "managed")
  job <- rrmg_output_begin(managed, "test")
  writeLines("old page", file.path(job$stage, "page_02.txt"))
  rrmg_output_finish(job)
  # A failed rendering invocation must leave the last completed output intact.
  job <- rrmg_output_begin(managed, "test")
  writeLines("partial", file.path(job$stage, "page_01.txt"))
  rrmg_output_abort(job)
  stopifnot(file.exists(file.path(managed, "page_02.txt")), !file.exists(file.path(managed, "page_01.txt")))
  job <- rrmg_output_begin(managed, "test")
  writeLines("new page", file.path(job$stage, "page_01.txt"))
  rrmg_output_finish(job)
  stopifnot(!file.exists(file.path(managed, "page_02.txt")))
  writeLines("edited", file.path(managed, "page_01.txt"))
  stopifnot(inherits(try(rrmg_output_begin(managed, "test"), silent = TRUE), "try-error"))
  stopifnot(readLines(file.path(managed, "page_01.txt")) == "edited")
  cat("[PASS] Failed render keeps previous output; stale pages removed; modified files protected.\n")
}
run_tests()
