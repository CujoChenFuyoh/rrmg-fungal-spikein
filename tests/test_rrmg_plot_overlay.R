#!/usr/bin/env Rscript
# Plot-layer regression with artificial profiles; never interpreted as study data.
# Usage: Rscript tests/test_rrmg_plot_overlay.R ../rrmg_plot_alignment.R OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
exprs <- parse(args[1])
is_assignment <- function(x) is.call(x) && identical(x[[1]], as.name("<-"))
main_expr <- Filter(function(x) is_assignment(x) && identical(x[[2]], as.name("main")), as.list(exprs))[[1]]
definitions <- as.list(main_expr[[3]][[3]])[-1]
wanted <- c("TARGET_READ_LABEL", "METAGENOME_READ_LABEL", "READ_SET_LEVELS",
            "READ_SET_COLORS", "READ_SET_FILL_COLORS", "DEPTH_TARGET_ALPHA", "DEPTH_SPIKED_ALPHA",
            "READ_BLOCK_ALPHA", "theme_gene", "add_fungamr_guides", "make_depth_panel", "make_read_blocks_panel")
for (definition in definitions) {
  if (is_assignment(definition) && is.symbol(definition[[2]]) && as.character(definition[[2]]) %in% wanted) eval(definition)
}
stopifnot(all(vapply(wanted, exists, logical(1))))
depth_fixture <- function(target, mixed) {
  rbind(data.frame(pos = seq_along(target), depth = target, ReadSource = TARGET_READ_LABEL),
        data.frame(pos = seq_along(mixed), depth = mixed, ReadSource = METAGENOME_READ_LABEL))
}
grey <- c(rep(0, 10), seq(0, 10, length.out = 20), rep(10, 20), seq(10, 0, length.out = 20), rep(0, 30))
fixtures <- list(identical = depth_fixture(grey, grey),
                 overlapping = depth_fixture(grey, c(rep(0, 30), grey[1:70] * 0.7)),
                 zero_mixed = depth_fixture(grey, rep(0, 100)))
plots <- lapply(fixtures, make_depth_panel, gene_len = 100, guide_df = data.frame())
for (p in plots) {
  stopifnot(length(p$layers) == 4L,
            p$layers[[1]]$aes_params$alpha == 0.20,
            p$layers[[2]]$aes_params$alpha == 0.45,
            p$layers[[3]]$aes_params$linetype == "dashed",
            p$layers[[4]]$aes_params$linetype == "solid",
            all(p$layers[[4]]$data$ReadSource == METAGENOME_READ_LABEL))
}
stopifnot(all(ggplot_build(plots$zero_mixed)$data[[4]]$y == 0),
          grepl("metagenomic reads: 0/100", plots$zero_mixed$labels$subtitle, fixed = TRUE))
blocks <- data.frame(target_start_1based = c(10, 10), target_end = c(80, 80), read_y = c(1, 5),
                     ReadSource = c(TARGET_READ_LABEL, METAGENOME_READ_LABEL))
block_plot <- make_read_blocks_panel(blocks, 100, guide_df = data.frame())
stopifnot(identical(ggplot_build(block_plot)$data[[1]]$alpha, c(0.5, 0.9)))
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
combined <- (plots$identical + labs(title = "Artificial test: identical profiles")) /
            (plots$overlapping + labs(title = "Artificial test: partial overlap")) /
            (plots$zero_mixed + labs(title = "Artificial test: zero mixed-read depth"))
ggsave(file.path(args[2], "depth_overlay_cases.png"), combined, width = 10.5, height = 9, dpi = 180)
ggsave(file.path(args[2], "depth_overlay_cases.pdf"), combined, width = 10.5, height = 9)
cat("[PASS] Explicit alpha/layer order, identical and partial overlap, zero preserved, separate read-track alpha.\n")
