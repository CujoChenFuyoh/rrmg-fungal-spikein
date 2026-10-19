#!/usr/bin/env Rscript
# Gene-alignment and reference-site coverage figures; native R or Windows R in WSL.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript rrmg_plot_alignment.R INPUT_DIR OUTPUT_DIR")
input_dir <- normalizePath(args[1], mustWork = TRUE)
if (!file.exists(file.path(input_dir, "COMPLETE.json"))) stop("Alignment inputs are incomplete; complete runbook section 14 first.")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file, mustWork = TRUE)), "rrmg_plot_common.R"))
# Wrap execution so errors clean up only this invocation's staging directory.
main <- function() {
job <- rrmg_output_begin(args[2], "rrmg_alignment_figures")
on.exit(rrmg_output_abort(job), add = TRUE)
base_out <- job$stage
# Labels are reference annotations/coverage, not variant or phenotype calls.


#===============================================================================
# CYP51A alignment and reference-site coverage from RRMG figure_data
#
# Read sets: mixed metagenome QC = orange; target-only synthetic reads = grey.
# Depth panels show overlapping coverage areas. Read blocks are exported as
# separate per-run and multi-page figures with a TSV alignment summary.
#
# Required input:
#   figure_data/CYP51A_alignment/plot_manifest.with_target_only.tsv
#   or figure_data/CYP51A_alignment/plot_manifest.tsv
#
# Mixed-read and contig inputs:
#   <PlotID>/contigs_to_GENE.paf
#   <PlotID>/reads_to_GENE.depth.tsv
#   <PlotID>/reads_to_GENE.paf.gz
#
# Target-only inputs:
#   <PlotID>/target_only_reads_to_GENE.depth.tsv
#   <PlotID>/target_only_reads_to_GENE.paf.gz
#
# Optional FungAMR mutation-site input:
#   figure_data/CYP51A_alignment/CYP51A_FungAMR_sites.tsv
#   Columns expected:
#     label, plot_start, plot_end, plot_mid, mutation_type, is_key_site
#===============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(tidyr)
library(purrr)

# Keep the existing ggplot/patchwork layout; export viewable PNG beside each PDF.
ggsave <- function(filename, plot, ...) {
  ggplot2::ggsave(filename, plot, ...)
  if (grepl("\\.pdf$", filename)) {
    ggplot2::ggsave(sub("\\.pdf$", ".png", filename), plot, ..., dpi = 180)
  }
}

#-------------------------------------------------------------------------------
# paths
#-------------------------------------------------------------------------------

figure_data_dir <- input_dir

manifest_file_preferred <- file.path(figure_data_dir, "plot_manifest.with_target_only.tsv")
manifest_file_default <- file.path(figure_data_dir, "plot_manifest.tsv")

manifest_file <- if (file.exists(manifest_file_preferred)) {
  manifest_file_preferred
} else {
  manifest_file_default
}

mutation_site_file <- file.path(figure_data_dir, "CYP51A_FungAMR_sites.tsv")

outdir <- file.path(base_out, "alignment_depth")
readblock_outdir <- file.path(base_out, "read_blocks")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(readblock_outdir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# settings
#-------------------------------------------------------------------------------

# Expected coverage is read from the manifest, not a hard-coded genome length.

MAX_READS_PER_SOURCE <- 350
MIN_READ_ALIGNMENT_LENGTH <- 20

# FungAMR mutation site display.
# "key": show all mutation ticks but label only key sites.
# "all": label all sites, may be crowded.
# "none": show ticks only.
FUNGAMR_LABEL_MODE <- "key"
FUNGAMR_KEY_REGEX <- "G54|L98H|Y121F|G138|M220|T289|G448S|Y491H|F495I"
FUNGAMR_SITE_COLOR <- "#CC79A7"
FUNGAMR_LABEL_SIZE <- 2.55

# Dashed vertical guide lines for annotated FungAMR sites across lower panels.
# "key": draw guides for key/labeled sites only; "all": draw every site; "none": no guides.
FUNGAMR_GUIDE_MODE <- "key"
FUNGAMR_GUIDE_COLOR <- "#B04A9F"
FUNGAMR_GUIDE_ALPHA <- 0.38
FUNGAMR_GUIDE_LINEWIDTH <- 0.23

# Read-set colors:
#   Spiked metagenome reads = orange, more visually prominent
#   Synthetic raw genome reads = grey baseline
READ_SET_LEVELS <- c("Synthetic raw genome reads", "Spiked metagenome reads")

READ_SET_COLORS <- c(
  "Synthetic raw genome reads" = "grey55",
  "Spiked metagenome reads" = "#D55E00"
)

READ_SET_FILL_COLORS <- c(
  "Synthetic raw genome reads" = "grey75",
  "Spiked metagenome reads" = "#D55E00"
)

# Alpha is opacity (0 = invisible, 1 = opaque). Keep mixed QC reads on top.
DEPTH_TARGET_ALPHA <- 0.20
DEPTH_SPIKED_ALPHA <- 0.45
READ_BLOCK_ALPHA <- c("Synthetic raw genome reads" = 0.50, "Spiked metagenome reads" = 0.90)

#-------------------------------------------------------------------------------
# helper functions
#-------------------------------------------------------------------------------

read_fasta_length <- function(fasta) {
  x <- readLines(fasta, warn = FALSE)
  sum(nchar(x[!grepl("^>", x)]))
}

short_label <- function(x, n = 42) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

wrap_text <- function(x, width = 90) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

safe_num <- function(x, default = NA_real_) {
  y <- suppressWarnings(as.numeric(x))
  ifelse(is.na(y), default, y)
}

load_fungamr_sites <- function(site_file) {
  if (!file.exists(site_file) || file.info(site_file, extra_cols = FALSE)$size == 0) {
    message("No FungAMR mutation-site file found: ", site_file)
    return(data.frame())
  }

  x <- read.delim(site_file, check.names = FALSE, stringsAsFactors = FALSE)

  required <- c("label", "plot_start", "plot_end")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop("Mutation-site file is missing required columns: ", paste(missing, collapse = ", "))
  }

  if (!"plot_mid" %in% colnames(x)) {
    x$plot_mid <- (as.numeric(x$plot_start) + as.numeric(x$plot_end)) / 2
  }

  if (!"is_key_site" %in% colnames(x)) {
    x$is_key_site <- ifelse(grepl(FUNGAMR_KEY_REGEX, x$label), "Y", "N")
  }

  if (!"mutation_type" %in% colnames(x)) {
    x$mutation_type <- "mutation"
  }

  x %>%
    mutate(
      plot_start = suppressWarnings(as.numeric(plot_start)),
      plot_end = suppressWarnings(as.numeric(plot_end)),
      plot_mid = suppressWarnings(as.numeric(plot_mid)),
      is_key_site = toupper(as.character(is_key_site)) %in% c("Y", "YES", "TRUE", "1"),
      label = as.character(label),
      mutation_type = as.character(mutation_type),
      label_root = ifelse(
        grepl("^[A-Za-z*][0-9]+", label),
        sub("^([A-Za-z*][0-9]+).*$", "\\1", label),
        label
      )
    ) %>%
    filter(!is.na(plot_mid), !is.na(plot_start), !is.na(plot_end))
}

make_fungamr_label_df <- function(site_df) {
  if (nrow(site_df) == 0 || FUNGAMR_LABEL_MODE == "none") {
    return(data.frame())
  }

  label_df <- site_df

  if (FUNGAMR_LABEL_MODE == "key") {
    label_df <- label_df %>% filter(is_key_site | grepl(FUNGAMR_KEY_REGEX, label))
  } else if (FUNGAMR_LABEL_MODE != "all") {
    stop("FUNGAMR_LABEL_MODE must be one of: key, all, none")
  }

  if (nrow(label_df) == 0) {
    return(data.frame())
  }

  # Collapse multiple labels at the same coordinate, e.g. G54E/K/R/V/W -> G54,
  # M220I/K/L/R/T/V/W -> M220. Singletons keep the full label.
  label_df %>%
    mutate(coord_key = round(plot_mid, 3)) %>%
    group_by(coord_key) %>%
    summarise(
      plot_mid = first(plot_mid),
      label = ifelse(n() > 1, first(label_root), first(label)),
      n_site = n(),
      .groups = "drop"
    ) %>%
    arrange(plot_mid) %>%
    mutate(
      label_track = (row_number() - 1) %% 3,
      label_y = 1.72 + 0.18 * label_track
    )
}


get_visible_fungamr_sites <- function(gene_len) {
  if (!exists("fungamr_sites") || nrow(fungamr_sites) == 0) {
    return(data.frame())
  }

  fungamr_sites %>%
    filter(
      !is.na(plot_mid),
      plot_mid >= 1,
      plot_mid <= gene_len,
      plot_end >= 1,
      plot_start <= gene_len
    ) %>%
    mutate(
      plot_start = pmax(1, as.numeric(plot_start)),
      plot_end = pmin(gene_len, as.numeric(plot_end)),
      plot_mid = as.numeric(plot_mid)
    )
}

make_fungamr_guide_df <- function(gene_len) {
  site_df <- get_visible_fungamr_sites(gene_len)

  if (nrow(site_df) == 0 || FUNGAMR_GUIDE_MODE == "none") {
    return(data.frame())
  }

  if (FUNGAMR_GUIDE_MODE == "key") {
    site_df <- site_df %>% filter(is_key_site | grepl(FUNGAMR_KEY_REGEX, label))
  } else if (FUNGAMR_GUIDE_MODE == "all") {
    site_df <- site_df
  } else {
    stop("FUNGAMR_GUIDE_MODE must be one of: key, all, none")
  }

  if (nrow(site_df) == 0) {
    return(data.frame())
  }

  site_df %>%
    mutate(coord_key = round(plot_mid, 3)) %>%
    group_by(coord_key) %>%
    summarise(
      plot_mid = first(plot_mid),
      label = ifelse(n() > 1, first(label_root), first(label)),
      n_site = n(),
      .groups = "drop"
    ) %>%
    arrange(plot_mid)
}

add_fungamr_guides <- function(p, guide_df) {
  if (nrow(guide_df) == 0) return(p)

  p +
    geom_vline(
      data = guide_df,
      aes(xintercept = plot_mid),
      linetype = "dashed",
      color = FUNGAMR_GUIDE_COLOR,
      linewidth = FUNGAMR_GUIDE_LINEWIDTH,
      alpha = FUNGAMR_GUIDE_ALPHA
    )
}

fungamr_sites <- load_fungamr_sites(mutation_site_file)
if (nrow(fungamr_sites) > 0) {
  message("Loaded FungAMR mutation sites: ", nrow(fungamr_sites))
}

read_paf <- function(paf_file) {
  if (is.na(paf_file) || !nzchar(paf_file)) {
    return(data.frame())
  }

  if (!file.exists(paf_file) || file.info(paf_file, extra_cols = FALSE)$size == 0) {
    return(data.frame())
  }

  con <- if (grepl("\\.gz$", paf_file)) {
    gzfile(paf_file, open = "rt")
  } else {
    file(paf_file, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  paf <- tryCatch(
    read.delim(
      con,
      header = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "",
      fill = TRUE
    ),
    error = function(e) data.frame()
  )

  if (nrow(paf) == 0 || ncol(paf) < 12) {
    return(data.frame())
  }

  paf <- paf[, 1:12]
  colnames(paf) <- c(
    "query_name", "query_len", "query_start", "query_end",
    "strand", "target_name", "target_len", "target_start",
    "target_end", "num_match", "aln_block_len", "mapq"
  )

  paf %>%
    mutate(
      target_start_1based = as.numeric(target_start) + 1,
      target_end = as.numeric(target_end),
      query_len = as.numeric(query_len),
      query_start = as.numeric(query_start),
      query_end = as.numeric(query_end),
      aln_block_len = as.numeric(aln_block_len),
      num_match = as.numeric(num_match),
      mapq = as.numeric(mapq),
      identity = ifelse(aln_block_len > 0, num_match / aln_block_len * 100, NA_real_),
      contig_label = query_name
    )
}

read_depth <- function(depth_file) {
  if (is.na(depth_file) || !nzchar(depth_file)) {
    return(data.frame())
  }

  if (!file.exists(depth_file) || file.info(depth_file, extra_cols = FALSE)$size == 0) {
    return(data.frame())
  }

  x <- read.delim(
    depth_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(x) == 0 || ncol(x) < 3) {
    return(data.frame())
  }

  colnames(x)[1:3] <- c("target_name", "pos", "depth")

  x %>%
    mutate(
      pos = as.numeric(pos),
      depth = as.numeric(depth)
    )
}

layout_read_tracks_one_source <- function(df, max_reads = 350, gap = 1) {
  if (nrow(df) == 0) return(df)

  df <- df %>%
    filter(!is.na(target_start_1based), !is.na(target_end)) %>%
    filter(aln_block_len >= MIN_READ_ALIGNMENT_LENGTH) %>%
    arrange(target_start_1based, target_end)

  if (nrow(df) == 0) return(df)

  df$n_total_source_reads <- nrow(df)

  if (nrow(df) > max_reads) {
    keep <- unique(round(seq(1, nrow(df), length.out = max_reads)))
    df <- df[keep, , drop = FALSE]
  }

  track_end <- numeric(0)
  track_id <- integer(nrow(df))

  for (i in seq_len(nrow(df))) {
    s <- df$target_start_1based[i]
    e <- df$target_end[i]
    placed <- FALSE

    if (length(track_end) > 0) {
      for (k in seq_along(track_end)) {
        if (s > track_end[k] + gap) {
          track_id[i] <- k
          track_end[k] <- e
          placed <- TRUE
          break
        }
      }
    }

    if (!placed) {
      track_end <- c(track_end, e)
      track_id[i] <- length(track_end)
    }
  }

  df$track_in_source <- track_id
  df$n_plotted_source_reads <- nrow(df)
  df
}

layout_combined_read_tracks <- function(spiked_paf, target_paf) {
  out <- list()

  if (nrow(target_paf) > 0) {
    x <- layout_read_tracks_one_source(target_paf, MAX_READS_PER_SOURCE)
    if (nrow(x) > 0) {
      x$ReadSource <- "Synthetic raw genome reads"
      out[["target_only"]] <- x
    }
  }

  if (nrow(spiked_paf) > 0) {
    x <- layout_read_tracks_one_source(spiked_paf, MAX_READS_PER_SOURCE)
    if (nrow(x) > 0) {
      x$ReadSource <- "Spiked metagenome reads"
      out[["spiked"]] <- x
    }
  }

  if (length(out) == 0) {
    return(data.frame())
  }

  df <- bind_rows(out)
  df$ReadSource <- factor(df$ReadSource, levels = READ_SET_LEVELS)

  max_target <- df %>%
    filter(ReadSource == "Synthetic raw genome reads") %>%
    summarise(m = ifelse(n() == 0, 0, max(track_in_source))) %>%
    pull(m)

  gap <- 3

  df <- df %>%
    mutate(
      read_y = ifelse(
        ReadSource == "Synthetic raw genome reads",
        track_in_source,
        track_in_source + max_target + gap
      )
    )

  df
}

summarise_read_blocks <- function(df, gene_len) {
  if (nrow(df) == 0) {
    return(data.frame(
      ReadSource = character(),
      NTotalAlignmentsBeforeSubsample = integer(),
      NPlottedAlignments = integer(),
      CoveredBpByBlocks = integer(),
      BreadthPctByBlocks = numeric(),
      MeanAlnBlockLen = numeric(),
      MedianAlnBlockLen = numeric(),
      MeanIdentity = numeric()
    ))
  }

  df %>%
    group_by(ReadSource) %>%
    summarise(
      NTotalAlignmentsBeforeSubsample = max(n_total_source_reads, na.rm = TRUE),
      NPlottedAlignments = n(),
      CoveredBpByBlocks = {
        intervals <- tibble(
          start = pmax(1, as.integer(target_start_1based)),
          end = pmin(gene_len, as.integer(target_end))
        ) %>%
          filter(!is.na(start), !is.na(end), end >= start) %>%
          arrange(start, end)

        if (nrow(intervals) == 0) {
          0L
        } else {
          merged_start <- intervals$start[1]
          merged_end <- intervals$end[1]
          total <- 0L

          if (nrow(intervals) > 1) {
            for (i in 2:nrow(intervals)) {
              if (intervals$start[i] <= merged_end + 1) {
                merged_end <- max(merged_end, intervals$end[i])
              } else {
                total <- total + merged_end - merged_start + 1
                merged_start <- intervals$start[i]
                merged_end <- intervals$end[i]
              }
            }
          }

          total <- total + merged_end - merged_start + 1
          total
        }
      },
      BreadthPctByBlocks = CoveredBpByBlocks / gene_len * 100,
      MeanAlnBlockLen = mean(aln_block_len, na.rm = TRUE),
      MedianAlnBlockLen = median(aln_block_len, na.rm = TRUE),
      MeanIdentity = mean(identity, na.rm = TRUE),
      .groups = "drop"
    )
}

theme_gene <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey40"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      plot.title = element_text(face = "bold", hjust = 0, lineheight = 0.95),
      plot.subtitle = element_text(hjust = 0, lineheight = 0.95),
      legend.position = "right"
    )
}

#-------------------------------------------------------------------------------
# data extraction for one run
#-------------------------------------------------------------------------------

get_run_plot_data <- function(row, base_dir) {
  target_fa <- file.path(base_dir, row$TargetFA)
  contig_paf <- file.path(base_dir, row$ContigPAF)
  depth_file <- file.path(base_dir, row$ReadDepthTSV)

  spiked_read_paf <- ""
  if ("ReadPAFGZ" %in% names(row) && !is.na(row$ReadPAFGZ) && nzchar(row$ReadPAFGZ)) {
    spiked_read_paf <- file.path(base_dir, row$ReadPAFGZ)
  }

  target_only_read_paf <- ""
  if ("TargetOnlyReadPAFGZ" %in% names(row) &&
      !is.na(row$TargetOnlyReadPAFGZ) &&
      nzchar(row$TargetOnlyReadPAFGZ)) {
    target_only_read_paf <- file.path(base_dir, row$TargetOnlyReadPAFGZ)
  }

  target_only_depth_file <- ""
  if ("TargetOnlyReadDepthTSV" %in% names(row) &&
      !is.na(row$TargetOnlyReadDepthTSV) &&
      nzchar(row$TargetOnlyReadDepthTSV)) {
    target_only_depth_file <- file.path(base_dir, row$TargetOnlyReadDepthTSV)
  }

  gene_len <- read_fasta_length(target_fa)

  paf <- read_paf(contig_paf)
  depth <- read_depth(depth_file)
  target_only_depth <- read_depth(target_only_depth_file)

  spiked_reads <- read_paf(spiked_read_paf)
  target_only_reads <- read_paf(target_only_read_paf)
  combined_reads <- layout_combined_read_tracks(spiked_reads, target_only_reads)

  if (nrow(depth) != gene_len || nrow(target_only_depth) != gene_len) {
    stop("Missing/incomplete depth data: ", row$RunID, "; missing is not zero.")
  }

  depth$ReadSource <- "Spiked metagenome reads"

  if (nrow(target_only_depth) > 0) {
    target_only_depth$ReadSource <- "Synthetic raw genome reads"
    depth_plot_df <- bind_rows(depth, target_only_depth)
  } else {
    depth_plot_df <- depth
  }

  depth_plot_df$ReadSource <- factor(depth_plot_df$ReadSource, levels = READ_SET_LEVELS)

  if (nrow(paf) > 0) {
    paf <- paf %>%
      arrange(target_start_1based, desc(aln_block_len)) %>%
      mutate(
        contig_label = factor(contig_label, levels = rev(unique(contig_label)))
      )
  }

  list(
    gene_len = gene_len,
    paf = paf,
    depth = depth_plot_df,
    depth_spiked = depth,
    depth_synthetic = target_only_depth,
    spiked_reads_raw = spiked_reads,
    synthetic_reads_raw = target_only_reads,
    combined_reads = combined_reads
  )
}

#-------------------------------------------------------------------------------
# plot panels
#-------------------------------------------------------------------------------

make_gene_panel <- function(row, gene_len) {
  run_id_line <- wrap_text(row$RunID, width = 70)

  info_line1 <- paste0(
    "ta=", row$TargetAbundanceDesign,
    ", total pairs=", format(as.numeric(row$TotalPairs), scientific = FALSE, big.mark = ","),
    ", seed=", row$Seed
  )

  target_only_breadth <- ""
  if ("TargetOnlyReadBreadthPct" %in% names(row) &&
      !is.na(row$TargetOnlyReadBreadthPct) &&
      nzchar(as.character(row$TargetOnlyReadBreadthPct))) {
    target_only_breadth <- paste0(
      ", synthetic read breadth=",
      round(as.numeric(row$TargetOnlyReadBreadthPct), 1),
      "%"
    )
  }

  info_line2 <- paste0(
    "Afum genome cov=", round(as.numeric(row$AfumGenomeCoverage), 3), "x",
    ", contig breadth=", round(as.numeric(row$ContigTargetBreadthPct), 1), "%",
    ", spiked read breadth=", round(as.numeric(row$ReadBreadthPct), 1), "%",
    target_only_breadth
  )

  title_text <- paste(run_id_line, info_line1, info_line2, sep = "\n")

  gene_df <- data.frame(
    xmin = 1,
    xmax = gene_len,
    ymin = 0.85,
    ymax = 1.15
  )

  site_df <- data.frame()
  label_df <- data.frame()

  if (exists("fungamr_sites") && nrow(fungamr_sites) > 0) {
    site_df <- fungamr_sites %>%
      filter(
        !is.na(plot_mid),
        plot_mid >= 1,
        plot_mid <= gene_len,
        plot_end >= 1,
        plot_start <= gene_len
      ) %>%
      mutate(
        plot_start = pmax(1, plot_start),
        plot_end = pmin(gene_len, plot_end)
      )

    label_df <- make_fungamr_label_df(site_df)
  }

  p <- ggplot(gene_df) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey35",
      color = "black",
      linewidth = 0.3
    )

  if (nrow(site_df) > 0) {
    p <- p +
      geom_segment(
        data = site_df,
        aes(x = plot_mid, xend = plot_mid, y = 1.15, yend = 1.48),
        inherit.aes = FALSE,
        color = FUNGAMR_SITE_COLOR,
        linewidth = 0.25,
        alpha = 0.75
      ) +
      geom_point(
        data = site_df,
        aes(x = plot_mid, y = 1.48),
        inherit.aes = FALSE,
        color = FUNGAMR_SITE_COLOR,
        size = 1.05,
        alpha = 0.85
      )
  }

  if (nrow(label_df) > 0) {
    p <- p +
      geom_segment(
        data = label_df,
        aes(x = plot_mid, xend = plot_mid, y = 1.50, yend = label_y - 0.03),
        inherit.aes = FALSE,
        color = FUNGAMR_SITE_COLOR,
        linewidth = 0.2,
        alpha = 0.55
      ) +
      geom_text(
        data = label_df,
        aes(x = plot_mid, y = label_y, label = label),
        inherit.aes = FALSE,
        color = FUNGAMR_SITE_COLOR,
        size = FUNGAMR_LABEL_SIZE,
        angle = 45,
        hjust = 0,
        vjust = 0.5,
        check_overlap = TRUE
      )
  }

  p +
    annotate(
      "text",
      x = gene_len / 2,
      y = 1.38,
      label = paste0("CYP51A reference gene (", gene_len, " bp)"),
      size = 3.4
    ) +
    scale_x_continuous(limits = c(1, gene_len), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, 3.0), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = title_text) +
    theme_gene() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 9.5, face = "bold", lineheight = 1.0)
    )
}

make_contig_panel <- function(paf, gene_len, guide_df = make_fungamr_guide_df(gene_len)) {
  if (nrow(paf) > 0) {
    p <- ggplot(paf)
    p <- add_fungamr_guides(p, guide_df)

    p +
      geom_segment(
        aes(
          x = target_start_1based,
          xend = target_end,
          y = contig_label,
          yend = contig_label
        ),
        linewidth = 1.3,
        color = "#0072B2",
        lineend = "round"
      ) +
      geom_point(
        aes(x = target_start_1based, y = contig_label),
        size = 1.2,
        color = "black"
      ) +
      scale_x_continuous(limits = c(1, gene_len), expand = c(0, 0)) +
      scale_y_discrete(labels = function(x) short_label(x, n = 26)) +
      labs(
        x = NULL,
        y = "Assembly contigs",
        subtitle = "Assembly contigs aligned to CYP51A"
      ) +
      theme_gene() +
      theme(
        plot.margin = margin(6, 10, 6, 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
  } else {
    p <- ggplot()
    p <- add_fungamr_guides(p, guide_df)

    p +
      annotate("text", x = gene_len / 2, y = 1, label = "No contig alignment", size = 4) +
      scale_x_continuous(limits = c(1, gene_len), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 2)) +
      labs(x = NULL, y = "Assembly contigs") +
      theme_gene() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
      )
  }
}

make_depth_panel <- function(depth_plot_df, gene_len, guide_df = make_fungamr_guide_df(gene_len)) {
  # Explicit layer order: pale target background, orange foreground, then outlines.
  # A zero mixed-read profile stays zero; covered-base counts disambiguate absence.
  target_depth <- depth_plot_df %>% filter(ReadSource == "Synthetic raw genome reads")
  spiked_depth <- depth_plot_df %>% filter(ReadSource == "Spiked metagenome reads")
  p <- ggplot(depth_plot_df, aes(x = pos, y = depth, color = ReadSource, fill = ReadSource)) +
    geom_area(
      data = target_depth,
      alpha = DEPTH_TARGET_ALPHA,
      colour = NA,
      position = "identity"
    ) +
    geom_area(
      data = spiked_depth,
      alpha = DEPTH_SPIKED_ALPHA,
      colour = NA,
      position = "identity"
    )

  p <- add_fungamr_guides(p, guide_df)

  p +
    # A wider dashed underlay remains visible even when both profiles coincide.
    geom_line(data = target_depth, linewidth = 0.80, alpha = 0.75, linetype = "dashed") +
    geom_line(data = spiked_depth, linewidth = 0.40, alpha = 1, linetype = "solid") +
    scale_color_manual(values = READ_SET_COLORS, name = "Read set", drop = FALSE) +
    scale_fill_manual(values = READ_SET_FILL_COLORS, name = "Read set", drop = FALSE) +
    scale_x_continuous(limits = c(1, gene_len), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      x = "Position on CYP51A reference (bp)",
      y = "Read depth",
      subtitle = sprintf("CYP51A covered bp | mixed QC: %d/%d; target-only: %d/%d",
                         sum(spiked_depth$depth > 0), gene_len, sum(target_depth$depth > 0), gene_len)
    ) +
    theme_gene() +
    theme(
      axis.text.x = element_text(),
      axis.ticks.x = element_line()
    )
}

make_read_blocks_panel <- function(combined_reads, gene_len, show_title = TRUE, guide_df = make_fungamr_guide_df(gene_len)) {
  if (nrow(combined_reads) > 0) {
    p <- ggplot(combined_reads)
    p <- add_fungamr_guides(p, guide_df)

    p +
      geom_segment(
        aes(
          x = target_start_1based,
          xend = target_end,
          y = read_y,
          yend = read_y,
          color = ReadSource,
          alpha = ReadSource
        ),
        linewidth = 0.55,
        lineend = "round"
      ) +
      scale_alpha_manual(values = READ_BLOCK_ALPHA, guide = "none") +
      scale_color_manual(values = READ_SET_COLORS, name = "Read set", drop = FALSE) +
      scale_x_continuous(limits = c(1, gene_len), expand = c(0, 0)) +
      scale_y_reverse(expand = expansion(mult = c(0.02, 0.06))) +
      labs(
        x = "Position on CYP51A (bp)",
        y = "Read blocks",
        subtitle = if (show_title) sprintf("Shown blocks | target-only: %d; mixed QC: %d (separate tracks)",
          sum(combined_reads$ReadSource == "Synthetic raw genome reads"),
          sum(combined_reads$ReadSource == "Spiked metagenome reads")) else NULL
      ) +
      theme_gene() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  } else {
    p <- ggplot()
    p <- add_fungamr_guides(p, guide_df)

    p +
      annotate("text", x = gene_len / 2, y = 1, label = "No read-level PAF available", size = 4) +
      scale_x_continuous(limits = c(1, gene_len), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 2)) +
      labs(x = "Position on CYP51A (bp)", y = "Read blocks") +
      theme_gene() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
      )
  }
}

plot_one_run <- function(row, base_dir) {
  dat <- get_run_plot_data(row, base_dir)

  p_gene <- make_gene_panel(row, dat$gene_len)
  p_contig <- make_contig_panel(dat$paf, dat$gene_len)
  p_depth <- make_depth_panel(dat$depth, dat$gene_len)

  # Main mechanism figure only keeps:
  #   1. CYP51A reference gene
  #   2. Assembly contig alignment
  #   3. Read-depth area plot
  # Read blocks are exported separately by plot_read_blocks_only().
  p_gene / p_contig / p_depth +
    plot_layout(heights = c(1.35, 1.1, 1.55), guides = "collect") +
    plot_annotation(caption = "Reference-site annotation and alignment-span coverage; not a mutation call.") &
    theme(legend.box = "horizontal", legend.position = "bottom")
}

plot_read_blocks_only <- function(row, base_dir) {
  dat <- get_run_plot_data(row, base_dir)

  title_text <- paste0(
    row$RunID,
    "\n",
    "ta=", row$TargetAbundanceDesign,
    ", total pairs=", format(as.numeric(row$TotalPairs), scientific = FALSE, big.mark = ","),
    ", seed=", row$Seed,
    ", Afum genome cov=", round(as.numeric(row$AfumGenomeCoverage), 3), "x"
  )

  p <- make_read_blocks_panel(dat$combined_reads, dat$gene_len, show_title = FALSE) +
    labs(title = wrap_text(title_text, width = 95)) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0, lineheight = 1.0)
    )

  p
}


#-------------------------------------------------------------------------------
# FungAMR site coverage summary helpers
#-------------------------------------------------------------------------------

calc_depth_stats_for_site <- function(depth_df, site_start, site_end) {
  site_start <- as.integer(round(site_start))
  site_end <- as.integer(round(site_end))
  if (is.na(site_start) || is.na(site_end) || site_end < site_start) {
    return(data.frame(
      DepthAtSiteMid = NA_real_, MeanDepthAcrossSite = NA_real_, MaxDepthAcrossSite = NA_real_,
      MinDepthAcrossSite = NA_real_, CoveredBpAcrossSite = NA_integer_, AnyCovered = "N", AllCovered = "N"
    ))
  }

  all_pos <- data.frame(pos = seq(site_start, site_end))

  x <- all_pos %>% left_join(depth_df %>% select(pos, depth), by = "pos")
  if (anyNA(x$depth)) stop("Incomplete depth at annotation site; missing is not zero.")

  mid <- as.integer(round((site_start + site_end) / 2))
  mid_depth <- x$depth[which.min(abs(x$pos - mid))]

  data.frame(
    DepthAtSiteMid = mid_depth,
    MeanDepthAcrossSite = mean(x$depth),
    MaxDepthAcrossSite = max(x$depth),
    MinDepthAcrossSite = min(x$depth),
    CoveredBpAcrossSite = sum(x$depth > 0),
    AnyCovered = ifelse(any(x$depth > 0), "Y", "N"),
    AllCovered = ifelse(all(x$depth > 0), "Y", "N")
  )
}

calc_alignment_stats_for_site <- function(aln_df, site_start, site_end, gene_len) {
  site_start <- as.integer(round(site_start))
  site_end <- as.integer(round(site_end))
  mid <- as.integer(round((site_start + site_end) / 2))

  if (nrow(aln_df) == 0 || is.na(site_start) || is.na(site_end)) {
    return(data.frame(
      NAlignmentsOverlappingSite = 0L,
      NAlignmentsCoveringSiteMid = 0L,
      CoveredBpByAlignments = 0L,
      AnyCoveredByAlignments = "N",
      AllCoveredByAlignments = "N",
      BestIdentityOverSite = NA_real_
    ))
  }

  hit <- aln_df %>%
    filter(
      !is.na(target_start_1based),
      !is.na(target_end),
      target_start_1based <= site_end,
      target_end >= site_start
    )

  mid_hit <- aln_df %>%
    filter(
      !is.na(target_start_1based),
      !is.na(target_end),
      target_start_1based <= mid,
      target_end >= mid
    )

  all_pos <- seq(site_start, site_end)
  covered <- rep(FALSE, length(all_pos))

  if (nrow(hit) > 0) {
    for (j in seq_len(nrow(hit))) {
      s <- max(site_start, as.integer(round(hit$target_start_1based[j])))
      e <- min(site_end, as.integer(round(hit$target_end[j])))
      if (!is.na(s) && !is.na(e) && e >= s) {
        covered[all_pos >= s & all_pos <= e] <- TRUE
      }
    }
  }

  best_identity <- NA_real_
  if (nrow(hit) > 0) {
    best_identity <- suppressWarnings(max(hit$identity, na.rm = TRUE))
    if (!is.finite(best_identity)) best_identity <- NA_real_
  }

  data.frame(
    NAlignmentsOverlappingSite = nrow(hit),
    NAlignmentsCoveringSiteMid = nrow(mid_hit),
    CoveredBpByAlignments = sum(covered),
    AnyCoveredByAlignments = ifelse(any(covered), "Y", "N"),
    AllCoveredByAlignments = ifelse(all(covered), "Y", "N"),
    BestIdentityOverSite = best_identity
  )
}

build_fungamr_site_summary_for_run <- function(row, dat) {
  site_df <- get_visible_fungamr_sites(dat$gene_len)
  if (nrow(site_df) == 0) {
    return(data.frame())
  }

  purrr::map_dfr(seq_len(nrow(site_df)), function(i) {
    site <- site_df[i, ]
    s <- as.integer(round(site$plot_start))
    e <- as.integer(round(site$plot_end))
    m <- as.integer(round(site$plot_mid))

    contig <- calc_alignment_stats_for_site(dat$paf, s, e, dat$gene_len)
    spiked_depth <- calc_depth_stats_for_site(dat$depth_spiked, s, e)
    synthetic_depth <- calc_depth_stats_for_site(dat$depth_synthetic, s, e)
    spiked_read_blocks <- calc_alignment_stats_for_site(dat$spiked_reads_raw, s, e, dat$gene_len)
    synthetic_read_blocks <- calc_alignment_stats_for_site(dat$synthetic_reads_raw, s, e, dat$gene_len)

    data.frame(
      PlotID = as.character(row$PlotID),
      RunID = as.character(row$RunID),
      Seed = as.character(row$Seed),
      TargetAbundanceDesign = as.numeric(row$TargetAbundanceDesign),
      TotalPairs = as.numeric(row$TotalPairs),
      AfumGenomeCoverage = as.numeric(row$AfumGenomeCoverage),
      SiteLabel = as.character(site$label),
      SiteLabelRoot = as.character(site$label_root),
      MutationType = as.character(site$mutation_type),
      IsKeySite = ifelse(isTRUE(site$is_key_site), "Y", "N"),
      Chrom = ifelse("chrom" %in% names(site), as.character(site$chrom), NA_character_),
      ChrStart = ifelse("chr_start" %in% names(site), as.character(site$chr_start), NA_character_),
      ChrEnd = ifelse("chr_end" %in% names(site), as.character(site$chr_end), NA_character_),
      PlotStart = s,
      PlotEnd = e,
      PlotMid = m,
      SiteWidthBp = e - s + 1,

      ContigNOverlappingSite = contig$NAlignmentsOverlappingSite,
      ContigNAtSiteMid = contig$NAlignmentsCoveringSiteMid,
      ContigCoveredBpAcrossSite = contig$CoveredBpByAlignments,
      ContigAnyCovered = contig$AnyCoveredByAlignments,
      ContigAllCovered = contig$AllCoveredByAlignments,
      ContigBestIdentityOverSite = contig$BestIdentityOverSite,

      SpikedReadDepthAtSiteMid = spiked_depth$DepthAtSiteMid,
      SpikedMeanDepthAcrossSite = spiked_depth$MeanDepthAcrossSite,
      SpikedMaxDepthAcrossSite = spiked_depth$MaxDepthAcrossSite,
      SpikedMinDepthAcrossSite = spiked_depth$MinDepthAcrossSite,
      SpikedCoveredBpAcrossSite = spiked_depth$CoveredBpAcrossSite,
      SpikedAnyReadCovered = spiked_depth$AnyCovered,
      SpikedAllReadCovered = spiked_depth$AllCovered,
      SpikedReadBlockNOverlappingSite = spiked_read_blocks$NAlignmentsOverlappingSite,
      SpikedReadBlockNAtSiteMid = spiked_read_blocks$NAlignmentsCoveringSiteMid,

      SyntheticReadDepthAtSiteMid = synthetic_depth$DepthAtSiteMid,
      SyntheticMeanDepthAcrossSite = synthetic_depth$MeanDepthAcrossSite,
      SyntheticMaxDepthAcrossSite = synthetic_depth$MaxDepthAcrossSite,
      SyntheticMinDepthAcrossSite = synthetic_depth$MinDepthAcrossSite,
      SyntheticCoveredBpAcrossSite = synthetic_depth$CoveredBpAcrossSite,
      SyntheticAnyReadCovered = synthetic_depth$AnyCovered,
      SyntheticAllReadCovered = synthetic_depth$AllCovered,
      SyntheticReadBlockNOverlappingSite = synthetic_read_blocks$NAlignmentsOverlappingSite,
      SyntheticReadBlockNAtSiteMid = synthetic_read_blocks$NAlignmentsCoveringSiteMid,

      SpikedMinusSyntheticDepthAtSiteMid = spiked_depth$DepthAtSiteMid - synthetic_depth$DepthAtSiteMid,
      SpikedMinusSyntheticMeanDepthAcrossSite = spiked_depth$MeanDepthAcrossSite - synthetic_depth$MeanDepthAcrossSite
    )
  })
}

#-------------------------------------------------------------------------------
# read manifest and plot
#-------------------------------------------------------------------------------

manifest <- read.delim(
  manifest_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

manifest <- manifest %>%
  mutate(
    TargetAbundanceDesign = as.numeric(TargetAbundanceDesign),
    TotalPairs = as.numeric(TotalPairs),
    TargetGenomeCoverageExpected = suppressWarnings(as.numeric(TargetGenomeCoverageExpected)),
    AfumGenomeCoverage = TargetGenomeCoverageExpected,
    TargetOnlyReadBreadthPct = if ("TargetOnlyReadBreadthPct" %in% names(.)) {
      suppressWarnings(as.numeric(TargetOnlyReadBreadthPct))
    } else {
      NA_real_
    }
  )

manifest$PlotID <- factor(manifest$PlotID, levels = manifest$PlotID)

message("Using manifest: ", manifest_file)
message("Runs to plot: ", nrow(manifest))

#-------------------------------------------------------------------------------
# export full mechanism figures and read-block-only figures
#-------------------------------------------------------------------------------

read_block_summary_all <- list()
fungamr_site_summary_all <- list()

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]

  dat <- get_run_plot_data(row, figure_data_dir)

  # Summary of read-block PAFs
  s <- summarise_read_blocks(dat$combined_reads, dat$gene_len)
  if (nrow(s) > 0) {
    s <- s %>%
      mutate(
        PlotID = as.character(row$PlotID),
        RunID = row$RunID,
        Seed = row$Seed,
        TargetAbundanceDesign = row$TargetAbundanceDesign,
        TotalPairs = row$TotalPairs,
        AfumGenomeCoverage = row$AfumGenomeCoverage,
        .before = 1
      )
    read_block_summary_all[[length(read_block_summary_all) + 1]] <- s
  }

  # Summary of FungAMR mutation-site coverage by contigs and reads
  site_s <- build_fungamr_site_summary_for_run(row, dat)
  if (nrow(site_s) > 0) {
    fungamr_site_summary_all[[length(fungamr_site_summary_all) + 1]] <- site_s
  }

  # Full figure
  # Main figure does NOT include the read-block layer.
  # It only contains gene model + assembly contigs + read-depth area plot.
  p_full <- {
    p_gene <- make_gene_panel(row, dat$gene_len)
    p_contig <- make_contig_panel(dat$paf, dat$gene_len)
    p_depth <- make_depth_panel(dat$depth, dat$gene_len)

    p_gene / p_contig / p_depth +
      plot_layout(heights = c(1.35, 1.1, 1.55), guides = "collect") +
      plot_annotation(caption = "Reference-site annotation and alignment-span coverage; not a mutation call.") &
      theme(legend.box = "horizontal", legend.position = "bottom")
  }

  prefix <- file.path(outdir, paste0(row$PlotID, "_CYP51A_alignment_depth_with_target_only_reads"))

  ggsave(
    paste0(prefix, ".pdf"),
    p_full,
    width = 10.5,
    height = 7.4
  )

  # Read-block-only figure
  p_blocks <- plot_read_blocks_only(row, figure_data_dir)

  block_prefix <- file.path(readblock_outdir, paste0(row$PlotID, "_CYP51A_read_blocks_only"))

  ggsave(
    paste0(block_prefix, ".pdf"),
    p_blocks,
    width = 10.5,
    height = 3.4
  )
}

#-------------------------------------------------------------------------------
# read-block summary TSV
#-------------------------------------------------------------------------------

if (length(read_block_summary_all) > 0) {
  read_block_summary <- bind_rows(read_block_summary_all)

  write.table(
    read_block_summary,
    file.path(readblock_outdir, "CYP51A_read_block_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

#-------------------------------------------------------------------------------
# FungAMR site coverage summary TSVs
#-------------------------------------------------------------------------------

if (length(fungamr_site_summary_all) > 0) {
  fungamr_site_summary <- bind_rows(fungamr_site_summary_all)

  write.table(
    fungamr_site_summary,
    file.path(base_out, "CYP51A_FungAMR_site_coverage_by_run.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  fungamr_site_summary_by_design <- fungamr_site_summary %>%
    group_by(
      SiteLabel, SiteLabelRoot, MutationType, IsKeySite,
      PlotStart, PlotEnd, PlotMid,
      TargetAbundanceDesign, TotalPairs
    ) %>%
    summarise(
      NRuns = n_distinct(RunID),
      MeanAfumGenomeCoverage = mean(AfumGenomeCoverage, na.rm = TRUE),
      ContigAnyCoveredRatePct = mean(ContigAnyCovered == "Y", na.rm = TRUE) * 100,
      ContigAllCoveredRatePct = mean(ContigAllCovered == "Y", na.rm = TRUE) * 100,
      MeanContigNAtSiteMid = mean(ContigNAtSiteMid, na.rm = TRUE),
      SpikedAnyReadCoveredRatePct = mean(SpikedAnyReadCovered == "Y", na.rm = TRUE) * 100,
      SpikedAllReadCoveredRatePct = mean(SpikedAllReadCovered == "Y", na.rm = TRUE) * 100,
      MeanSpikedReadDepthAtSiteMid = mean(SpikedReadDepthAtSiteMid, na.rm = TRUE),
      MedianSpikedReadDepthAtSiteMid = median(SpikedReadDepthAtSiteMid, na.rm = TRUE),
      MeanSpikedMeanDepthAcrossSite = mean(SpikedMeanDepthAcrossSite, na.rm = TRUE),
      SyntheticAnyReadCoveredRatePct = mean(SyntheticAnyReadCovered == "Y", na.rm = TRUE) * 100,
      SyntheticAllReadCoveredRatePct = mean(SyntheticAllReadCovered == "Y", na.rm = TRUE) * 100,
      MeanSyntheticReadDepthAtSiteMid = mean(SyntheticReadDepthAtSiteMid, na.rm = TRUE),
      MedianSyntheticReadDepthAtSiteMid = median(SyntheticReadDepthAtSiteMid, na.rm = TRUE),
      MeanSyntheticMeanDepthAcrossSite = mean(SyntheticMeanDepthAcrossSite, na.rm = TRUE),
      MeanSpikedMinusSyntheticDepthAtSiteMid = mean(SpikedMinusSyntheticDepthAtSiteMid, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(TargetAbundanceDesign, TotalPairs, PlotMid, SiteLabel)

  write.table(
    fungamr_site_summary_by_design,
    file.path(base_out, "CYP51A_FungAMR_site_coverage_by_design.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

#-------------------------------------------------------------------------------
# multi-page PDFs
#-------------------------------------------------------------------------------

pdf(
  file.path(outdir, "CYP51A_alignment_depth_with_target_only_reads_all_runs.pdf"),
  width = 10.5,
  height = 7.4,
  onefile = TRUE
)

for (i in seq_len(nrow(manifest))) {
  print(plot_one_run(manifest[i, ], figure_data_dir))
}

dev.off()

pdf(
  file.path(readblock_outdir, "CYP51A_read_blocks_only_all_runs.pdf"),
  width = 10.5,
  height = 3.4,
  onefile = TRUE
)

for (i in seq_len(nrow(manifest))) {
  print(plot_read_blocks_only(manifest[i, ], figure_data_dir))
}

dev.off()

writeLines(c(
  "Reference site labels are annotations only, not detected mutations or AMR phenotypes.",
  "Site inclusion/exclusion and reference validation: see INPUT_DIR/site_annotation_audit.tsv.",
  "Main figures: gene model, assembly blocks, mixed/target-only depth.",
  "Read-block figures: a display-limited selection, not all input reads; see read-block summary.",
  "Coverage summaries describe aligned reference spans/depth, not base or codon correctness.",
  "Expected target genome coverage comes from the design, not actual FASTQ counts."
), file.path(base_out, "README.txt"))
capture.output(sessionInfo(), file = file.path(base_out, "R_sessionInfo.txt"))
rrmg_output_finish(job)
message("Done. Full figures: ", file.path(job$final, "alignment_depth"))
message("Read blocks: ", file.path(job$final, "read_blocks"))
message("Site coverage summaries (if annotations present): ", job$final)
}
main()
