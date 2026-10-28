#!/usr/bin/env Rscript
# RRMG local preview; base R only. No packages, reads regeneration or PBS needed.
# Usage: Rscript rrmg_plot_preview.R WORKDIR OUTPUT_DIR
# Load the bundled output helper from this script's directory, not the working directory.
.rrmg_preview_file <- if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(.rrmg_preview_file, mustWork = TRUE)), "rrmg_plot_common.R"))
# Read-only inputs: plan, three collected summary TSVs and per-run gene depths.
# Figure families follow the existing RRMG permutation-summary and target-genome
# plotting scripts, with portable paths and explicit missing-result denominators.
# Fractions are descriptive simulation output, not mutation calls or sensitivity estimates.

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             colClasses = "character", check.names = FALSE,
             na.strings = c("", "NA", "NaN", "nan"), stringsAsFactors = FALSE)
}

nonempty <- function(path) {
  file.exists(path) && !dir.exists(path) && file.info(path, extra_cols = FALSE)$size > 0
}

require_columns <- function(x, columns, name) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(name, ": missing columns: ", paste(missing, collapse = ", "))
}

check_ids <- function(x, field, name) {
  ids <- x[[field]]
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(name, ": empty or duplicate ", field, "; recollect/repair before plotting.")
  }
}

numeric_field <- function(x, field, name, maximum = Inf) {
  raw <- x[[field]]
  value <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & (!is.finite(value) | value < 0 | value > maximum)
  if (any(bad)) stop(name, ": invalid ", field, " (e.g. ", raw[which(bad)[1]], ")")
  value
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

collect_preview <- function(wd) {
  notes <- character()
  note <- function(x) notes <<- c(notes, x)
  plan_file <- file.path(wd, "input/rrmg/rrmg_permutation_plan.tsv")
  if (!nonempty(plan_file)) {
    plan_file <- file.path(wd, "figure_data/metadata/selected_plan.tsv")
    note("Using figure_data/metadata/selected_plan.tsv: figures describe this selected cohort only.")
  }
  plan <- read_tsv(plan_file)
  require_columns(plan, c("run_id", "run_index", "seed", "target_abundance", "total_pairs"), "plan")
  check_ids(plan, "run_id", "plan")
  check_ids(plan, "run_index", "plan")
  if (!nrow(plan)) stop("Plan has no samples.")
  if (any(grepl("[/\\\\]", plan$run_id) | plan$run_id %in% c(".", ".."))) {
    stop("Invalid run_id: expected a directory name, not a path.")
  }
  ids <- plan$run_id
  n <- length(ids)
  # Plan order is authoritative. Never join a result to a sample by row number.
  data <- data.frame(RunID = ids, RowIndex = plan$run_index, Seed = plan$seed,
                     TargetAbundanceDesign = numeric_field(plan, "target_abundance", "plan", 1),
                     TotalPairsDesign = numeric_field(plan, "total_pairs", "plan"),
                     stringsAsFactors = FALSE)
  if (anyNA(data[c("TargetAbundanceDesign", "TotalPairsDesign")])) {
    stop("Plan has missing abundance or total read-pair values; cannot group every planned sample.")
  }
  data$Label <- paste0("row ", data$RowIndex, " / seed ", data$Seed)
  status <- data[c("RunID", "RowIndex", "Seed")]
  specs <- list(
    contig = list(file = "rrmg_target_gene_map_merged.tsv", folder = "contig_mapping",
                  fields = c(ContigGeneBreadthPct = "TargetBreadthPct")),
    reads = list(file = "rrmg_reads_to_gene_depth_summary.tsv", folder = "read_to_gene",
                 fields = c(ReadGeneBreadthPct = "BreadthPct", ReadGeneMeanDepth = "MeanDepth",
                            ReadGeneBp = "TargetGeneBp")),
    target = list(file = "rrmg_afum_target_reads_genome_coverage_summary.tsv", folder = "target_genome",
                  fields = c(TargetGenomeBreadthPct = "ObservedBreadthPct",
                             TargetGenomeMeanDepth = "ObservedMeanDepth",
                             TargetPairsDesign = "TargetPairsDesign", TargetPairsFastq = "TargetPairsFastq"))
  )
  tables <- list()
  input_files <- plan_file
  for (kind in names(specs)) {
    spec <- specs[[kind]]
    path <- file.path(wd, "result", spec$file)
    if (!nonempty(path)) path <- file.path(wd, "figure_data", spec$folder, spec$file)
    table <- data.frame(RunID = character(), stringsAsFactors = FALSE)
    for (field in unname(spec$fields)) table[[field]] <- character()
    if (nonempty(path)) {
      table <- read_tsv(path)
      require_columns(table, c("RunID", unname(spec$fields)), spec$file)
      check_ids(table, "RunID", spec$file)
      input_files <- c(input_files, path)
      extras <- setdiff(table$RunID, ids)
      if (length(extras)) note(paste0(spec$file, ": ignored ", length(extras), " RunID(s) outside this plan."))
    } else {
      note(paste0(spec$file, ": missing/empty; values remain NA."))
    }
    index <- match(ids, table$RunID)
    status[[paste0(kind, "_summary")]] <- ifelse(is.na(index), "missing", "present")
    missing <- ids[is.na(index)]
    if (length(missing)) note(paste0(kind, " summaries missing for: ", paste(missing, collapse = ", ")))
    for (output in names(spec$fields)) {
      maximum <- if (grepl("Pct$", output)) 100 else Inf
      values <- numeric_field(table, spec$fields[[output]], spec$file, maximum)
      data[[output]] <- values[index]
      if (anyNA(data[[output]])) note(paste0(output, ": ", sum(is.na(data[[output]])), " missing value(s); not replaced with zero."))
    }
    tables[[kind]] <- table[index, , drop = FALSE]
  }
  # Additional metrics used in the original grouped plots. Absent columns stay NA.
  extra_fields <- list(contig = c(ContigCoveredBp = "TargetCoveredBp", ContigBestAlnLen = "BestAlnLen",
                                   ContigBestIdentityPct = "BestIdentity"),
                       reads = c(ReadMappedAlignments = "MappedReadAlignments", ReadMaxDepth = "MaxDepth",
                                  ReadBreadthGE5Pct = "BreadthPctDepthGE5"),
                       target = c(TargetCovFromFastq = "TargetCovFromFastq"))
  for (kind in names(extra_fields)) {
    for (column in names(extra_fields[[kind]])) {
      field <- extra_fields[[kind]][[column]]
      data[[column]] <- if (field %in% names(tables[[kind]]))
        numeric_field(tables[[kind]], field, kind, if (grepl("Pct$", column)) 100 else Inf) else rep(NA_real_, n)
    }
  }
  hit <- if ("HasGeneHit" %in% names(tables$contig)) toupper(tables$contig$HasGeneHit) else rep(NA_character_, n)
  if (any(!is.na(hit) & !hit %in% c("Y", "N"))) stop("HasGeneHit must be Y, N or NA.")
  data$ContigHit <- ifelse(is.na(hit), NA, hit == "Y")
  data$ReadHit <- data$ReadMappedAlignments > 0 | data$ReadGeneBreadthPct > 0
  data$BothHit <- ifelse(is.na(data$ContigHit) | is.na(data$ReadHit), NA, data$ContigHit & data$ReadHit)
  data$GenomeDepthObservedOverExpected <- ifelse(data$TargetCovFromFastq > 0,
                                                 data$TargetGenomeMeanDepth / data$TargetCovFromFastq, NA_real_)
  measured <- c("ContigGeneBreadthPct", "ReadGeneBreadthPct", "ReadGeneMeanDepth",
                "TargetGenomeBreadthPct", "TargetPairsFastq")
  if (all(is.na(as.matrix(data[measured])))) stop("No numeric results for this plan; collect results before plotting.")

  data$TargetPairsActualOverDesign <- ifelse(data$TargetPairsDesign > 0,
                                            data$TargetPairsFastq / data$TargetPairsDesign, NA_real_)
  known <- !is.na(data$TargetPairsDesign) & !is.na(data$TargetPairsFastq)
  status$target_pairs <- ifelse(!known, "missing",
                                ifelse(data$TargetPairsDesign == data$TargetPairsFastq, "match", "MISMATCH"))
  if (any(status$target_pairs == "MISMATCH")) {
    note(paste0("Target pairs differ from design: ", paste(ids[status$target_pairs == "MISMATCH"], collapse = ", "),
                ". Inspect read simulation; design abundance is not measured abundance."))
  }

  depth_parts <- list()
  panels <- list()
  status$gene_depth <- "missing"
  status$depth_vs_summary <- "not_checked"
  for (i in seq_len(n)) {
    # Prefer the relocated working tree over absolute paths saved on an old HPC.
    path <- file.path(wd, "temp/rrmg_permutation", ids[i], "readmap/reads_to_GENE.depth.tsv")
    if (!nonempty(path)) {
      portable <- file.path(
        wd, "figure_data/CYP51A_alignment",
        sprintf("row_%05d", as.integer(data$RowIndex[i])),
        "reads_to_GENE.depth.tsv"
      )
      if (nonempty(portable)) path <- portable
    }
    if (!nonempty(path) && "DepthTSV" %in% names(tables$reads)) {
      alternate <- tables$reads$DepthTSV[i]
      if (!is.na(alternate) && nzchar(alternate)) {
        if (!grepl("^(/|[A-Za-z]:[/\\\\])", alternate)) alternate <- file.path(wd, alternate)
        if (nonempty(alternate)) path <- alternate
      }
    }
    if (!nonempty(path)) {
      note(paste0(ids[i], ": gene depth file missing; no zero-depth curve inferred."))
      panels[[length(panels) + 1L]] <- list(i = i, ref = "missing depth", data = NULL)
      next
    }
    d <- tryCatch(read.delim(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                             colClasses = "character", stringsAsFactors = FALSE), error = function(e) NULL)
    valid <- !is.null(d) && ncol(d) == 3L && nrow(d) > 0L
    if (valid) {
      names(d) <- c("Reference", "Position", "Depth")
      d$Position <- suppressWarnings(as.numeric(d$Position))
      d$Depth <- suppressWarnings(as.numeric(d$Depth))
      valid <- !anyNA(d) && all(nzchar(d$Reference)) &&
        all(is.finite(d$Position) & d$Position >= 1 & d$Position == floor(d$Position)) &&
        all(is.finite(d$Depth) & d$Depth >= 0 & d$Depth == floor(d$Depth)) &&
        !anyDuplicated(d[c("Reference", "Position")])
    }
    if (!valid) {
      status$gene_depth[i] <- "invalid"
      note(paste0(ids[i], ": invalid depth file; expected unique reference/position/depth rows."))
      panels[[length(panels) + 1L]] <- list(i = i, ref = "invalid depth", data = NULL)
      next
    }
    status$gene_depth[i] <- "present"
    input_files <- c(input_files, path)
    pieces <- split(d, factor(d$Reference, levels = unique(d$Reference)))
    for (piece in pieces) {
      piece <- piece[order(piece$Position), , drop = FALSE]
      if (piece$Position[1] != 1 || any(diff(piece$Position) != 1)) {
        status$gene_depth[i] <- "gaps"
        note(paste0(ids[i], "/", piece$Reference[1], ": depth positions have gaps; absent bases are not zero-filled."))
      }
      panels[[length(panels) + 1L]] <- list(i = i, ref = piece$Reference[1], data = piece)
    }
    depth_parts[[length(depth_parts) + 1L]] <- data.frame(RunID = ids[i], d, SourceFile = path)
    if (!is.na(data$ReadGeneBp[i]) && data$ReadGeneBp[i] != nrow(d)) {
      status$depth_vs_summary[i] <- "MISMATCH"
    } else if (!is.na(data$ReadGeneMeanDepth[i]) && !is.na(data$ReadGeneBp[i])) {
      status$depth_vs_summary[i] <- ifelse(abs(mean(d$Depth) - data$ReadGeneMeanDepth[i]) < 1e-4, "match", "MISMATCH")
    }
    if (status$depth_vs_summary[i] == "MISMATCH") note(paste0(ids[i], ": depth rows/mean disagree with summary; check for stale or partial files."))
  }
  depth <- if (length(depth_parts)) do.call(rbind, depth_parts) else
    data.frame(RunID = character(), Reference = character(), Position = numeric(), Depth = numeric(), SourceFile = character())
  list(data = data, status = status, depth = depth, panels = panels,
       notes = unique(notes), inputs = unique(input_files), wd = wd)
}

# Paired point plots retain a true zero as a point, while NA is text, not a zero.
point_panel <- function(values, labels, title, ylabel, series, fixed_percent = FALSE) {
  values <- as.matrix(values)
  colors <- c("#287C8E", "#B16B31")
  marks <- c(16, 17)
  n <- nrow(values)
  k <- ncol(values)
  upper <- if (fixed_percent) 100 else max(c(0, values), na.rm = TRUE)
  if (upper == 0) upper <- 1
  plot(NA, xlim = c(0.4, n + 0.6), ylim = c(0, upper * 1.13), xaxs = "i", yaxs = "i",
       xaxt = "n", xlab = "", ylab = ylabel, main = title, bty = "l", cex.main = 0.95)
  abline(h = pretty(c(0, upper)), col = "#E7EAEC", lwd = 0.6)
  axis(1, at = seq_len(n), labels = labels, las = 2, cex.axis = 0.72, tick = FALSE)
  offsets <- if (k == 1) 0 else seq(-0.13, 0.13, length.out = k)
  for (j in seq_len(k)) {
    x <- seq_len(n) + offsets[j]
    points(x, values[, j], pch = marks[j], col = colors[j], cex = 1.05)
    missing <- is.na(values[, j])
    if (any(missing)) text(x[missing], rep(upper * (0.07 + 0.07 * (j - 1)), sum(missing)),
                           labels = "NA", col = colors[j], cex = 0.65)
  }
  legend("topleft", legend = series, col = colors[seq_len(k)], pch = marks[seq_len(k)],
         bty = "n", cex = 0.75, horiz = TRUE, inset = c(0, -0.01))
}

draw_overview <- function(bundle, index, page, total) {
  d <- bundle$data[index, , drop = FALSE]
  par(mfrow = c(2, 2), mar = c(6.5, 4.1, 2.4, 0.9), oma = c(3.3, 0.2, 3.7, 0.2),
      family = "sans", las = 1, mgp = c(2.5, 0.8, 0), cex = 0.9)
  point_panel(d[c("ContigGeneBreadthPct", "ReadGeneBreadthPct")], d$Label,
              "a  Target-gene breadth", "Covered reference bases (%)", c("Contigs", "Mixed QC reads"), TRUE)
  point_panel(d["ReadGeneMeanDepth"], d$Label, "b  Target-gene read depth", "Mean depth (x)", "Mixed QC reads")
  point_panel(d["TargetGenomeBreadthPct"], d$Label, "c  Target-only genome breadth",
              "Covered reference bases (%)", "Target-only reads", TRUE)
  point_panel(d[c("TargetPairsDesign", "TargetPairsFastq")] / 1000, d$Label,
              "d  Target-only read-pair audit", "Read pairs (thousands)", c("Designed", "FASTQ count"))
  mtext("RRMG | per-run output preview", side = 3, outer = TRUE, line = 2, font = 2, cex = 1.2)
  mtext(sprintf("%d planned simulation runs | page %d/%d | no significance tests", nrow(bundle$data), page, total),
        side = 3, outer = TRUE, line = 0.6, cex = 0.83)
  mismatches <- sum(bundle$status$target_pairs == "MISMATCH")
  mtext(sprintf("Target-pair mismatches: %d. NA = missing, not zero. See sample_status.tsv and warnings.txt.", mismatches),
        side = 1, outer = TRUE, line = 1, cex = 0.76)
  mtext("Workflow QC only; not a detection-rate estimate, mutation call, or measured species abundance.",
        side = 1, outer = TRUE, line = 2.2, cex = 0.76)
}

draw_depth <- function(bundle, panels, page, total) {
  par(mfrow = c(length(panels), 1), mar = c(3.8, 4.2, 2.6, 0.9),
      oma = c(2.4, 0.2, 3.5, 0.2), family = "sans", las = 1, mgp = c(2.1, 0.7, 0), cex = 0.9)
  for (panel in panels) {
    label <- paste0(bundle$data$Label[panel$i], " | ", panel$ref)
    if (is.null(panel$data)) {
      plot.new()
      title(main = label, cex.main = 0.95)
      text(0.5, 0.5, "Depth data missing or invalid (not zero)", col = "#777777")
      next
    }
    d <- panel$data
    plot(NA, xlim = c(1, max(2, d$Position)), ylim = c(0, max(1, d$Depth) * 1.08),
         xaxs = "i", yaxs = "i", bty = "l", xlab = "Reference position (bp)", ylab = "Read depth (x)",
         main = label, cex.main = 0.95)
    abline(h = pretty(c(0, max(1, d$Depth))), col = "#E7EAEC", lwd = 0.6)
    # Separate reference records and disconnected position segments; never bridge missing bases.
    segments <- split(seq_len(nrow(d)), cumsum(c(TRUE, diff(d$Position) != 1)))
    for (segment in segments) {
      if (length(segment) == 1) points(d$Position[segment], d$Depth[segment], pch = 16, cex = 0.5, col = "#287C8E")
      else lines(d$Position[segment], d$Depth[segment], col = "#287C8E", lwd = 1)
    }
  }
  mtext("RRMG | mixed QC reads to target gene", side = 3, outer = TRUE, line = 1.8, font = 2, cex = 1.2)
  mtext(sprintf("Per-base depth, unsmoothed | reference records shown separately | page %d/%d", page, total),
        side = 3, outer = TRUE, line = 0.5, cex = 0.83)
  mtext("Source: samtools depth TSVs. Missing positions are not filled with zero; axes may differ between panels.",
        side = 1, outer = TRUE, line = 0.9, cex = 0.76)
}

design_summary <- function(d) {
  groups <- split(seq_len(nrow(d)), interaction(d$TargetAbundanceDesign, d$TotalPairsDesign, drop = TRUE))
  answer <- lapply(groups, function(index) {
    g <- d[index, , drop = FALSE]
    out <- data.frame(TargetAbundanceDesign = g$TargetAbundanceDesign[1], TotalPairsDesign = g$TotalPairsDesign[1],
                      PlannedN = nrow(g), stringsAsFactors = FALSE)
    for (kind in c("Contig", "Read", "Both")) {
      hit <- g[[paste0(kind, "Hit")]]
      valid <- sum(!is.na(hit))
      out[[paste0(kind, "ValidN")]] <- valid
      out[[paste0(kind, "MissingN")]] <- sum(is.na(hit))
      out[[paste0(kind, "HitN")]] <- sum(hit, na.rm = TRUE)
      out[[paste0(kind, "HitPct")]] <- if (valid > 0) 100 * sum(hit, na.rm = TRUE) / valid else NA_real_
    }
    out
  })
  if (!length(answer)) stop("Design abundance/read-pair values are missing; cannot group results.")
  out <- do.call(rbind, answer)
  out[order(out$TargetAbundanceDesign, out$TotalPairsDesign), ]
}

pair_label <- function(x) {
  ifelse(x >= 1e6, paste0(format(x / 1e6, trim = TRUE, scientific = FALSE), "M"),
         ifelse(x >= 1e3, paste0(format(x / 1e3, trim = TRUE, scientific = FALSE), "K"), as.character(x)))
}

draw_hit_rates <- function(d) {
  abundance <- sort(unique(d$TargetAbundanceDesign))
  pairs <- sort(unique(d$TotalPairsDesign))
  par(mfrow = c(1, 2), mar = c(4.3, 4.8, 2.9, 1), oma = c(3.1, 0.2, 2.3, 0.2),
      family = "sans", las = 1, mgp = c(2.8, 0.7, 0))
  for (kind in c("Contig", "Read")) {
    palette <- colorRampPalette(c("#F1F3F4", if (kind == "Contig") "#0072B2" else "#D55E00"))(101)
    plot(NA, xlim = c(0.5, length(pairs) + 0.5), ylim = c(0.5, length(abundance) + 0.5),
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
         xlab = "Designed total metagenomic read pairs", ylab = "Designed target read fraction (%)",
         main = paste(kind, "target-gene hit fraction"), cex.main = 0.95)
    axis(1, at = seq_along(pairs), labels = pair_label(pairs), tick = FALSE)
    axis(2, at = seq_along(abundance), labels = format(100 * abundance, trim = TRUE), tick = FALSE)
    for (i in seq_len(nrow(d))) {
      x <- match(d$TotalPairsDesign[i], pairs)
      y <- match(d$TargetAbundanceDesign[i], abundance)
      rate <- d[[paste0(kind, "HitPct")]][i]
      col <- if (is.na(rate)) "#DDDDDD" else palette[1 + round(rate)]
      rect(x - 0.48, y - 0.48, x + 0.48, y + 0.48, col = col, border = "white")
      label <- sprintf("%d/%d valid\n%s | %d planned", d[[paste0(kind, "HitN")]][i],
                        d[[paste0(kind, "ValidN")]][i], if (is.na(rate)) "NA" else sprintf("%.1f%%", rate), d$PlannedN[i])
      text(x, y, label, cex = 0.76, col = if (!is.na(rate) && rate > 65) "white" else "#222222")
    }
  }
  mtext("RRMG | descriptive hit fractions by design", side = 3, outer = TRUE, line = 0.7, font = 2, cex = 1.1)
  mtext("Denominator = runs with evaluable hit status, not all planned runs; missing results excluded.",
        side = 1, outer = TRUE, line = 0.9, cex = 0.8)
  mtext("Contigs: upstream HasGeneHit. Reads: mapped alignments > 0 OR breadth > 0. No sensitivity claim.",
        side = 1, outer = TRUE, line = 2, cex = 0.8)
}

draw_grouped_metrics <- function(d, metrics, abundances, title) {
  pairs <- sort(unique(d$TotalPairsDesign))
  par(mfrow = c(length(metrics), length(abundances)), mar = c(3.8, 4.7, 2.2, 0.8),
      oma = c(3.2, 0, 3, 0), family = "sans", las = 1, mgp = c(2.0, 0.65, 0), cex = 0.85)
  for (field in names(metrics)) {
    finite <- d[[field]][is.finite(d[[field]])]
    upper <- if (grepl("Pct$", field)) 100 else max(c(1, finite))
    for (abundance in abundances) {
      g <- d[d$TargetAbundanceDesign == abundance, , drop = FALSE]
      values <- g[[field]]
      x <- match(g$TotalPairsDesign, pairs)
      plot(NA, xlim = c(0.45, length(pairs) + 0.55), ylim = c(0, upper * 1.18),
           xaxs = "i", yaxs = "i", xaxt = "n", bty = "l", xlab = "Total read pairs (design)",
           ylab = metrics[[field]], main = sprintf("Target read fraction: %g%%", 100 * abundance), cex.main = 0.95)
      abline(h = pretty(c(0, upper)), col = "#E7EAEC", lwd = 0.6)
      if (field == "GenomeDepthObservedOverExpected") abline(h = 1, col = "#777777", lty = 2)
      axis(1, at = seq_along(pairs), labels = pair_label(pairs), cex.axis = 0.8)
      for (j in seq_along(pairs)) {
        index <- which(x == j & is.finite(values))
        if (length(index) >= 3) {
          boxplot(values[index], at = j, add = TRUE, axes = FALSE, outline = FALSE,
                  boxwex = 0.6, col = "#EEEEEE", border = "#777777", lwd = 0.6)
        }
        if (length(index)) {
          # Deterministic horizontal spread; repeated plotting does not move points randomly.
          offset <- if (length(index) == 1) 0 else seq(-0.15, 0.15, length.out = length(index))
          color <- ifelse(is.na(g$ContigHit[index]), "#B16B31", ifelse(g$ContigHit[index], "#0072B2", "#888888"))
          points(j + offset, values[index], pch = 16, cex = 0.65, col = adjustcolor(color, alpha.f = 0.7))
        }
        text(j, upper * 1.1, sprintf("n=%d/%d", length(index), sum(x == j)), cex = 0.65, col = "#555555")
      }
    }
  }
  mtext(title, side = 3, outer = TRUE, line = 1.1, font = 2, cex = 1.15)
  mtext("Each dot = one simulation run. n = available/planned for each metric. Boxes only for n >= 3.",
        side = 1, outer = TRUE, line = 0.9, cex = 0.78)
  mtext("Contig hit: blue; no hit: grey; status missing: brown. Missing values are not zero-filled; no tests.",
        side = 1, outer = TRUE, line = 2, cex = 0.78)
}

draw_validation <- function(d) {
  par(mfrow = c(1, 2), mar = c(4.9, 4.8, 3, 1), oma = c(2.8, 0.2, 2.5, 0.2),
      family = "sans", las = 1, mgp = c(2.9, 0.8, 0), cex = 0.9)
  fields <- list(c("TargetCovFromFastq", "TargetGenomeMeanDepth"), c("TargetPairsDesign", "TargetPairsFastq"))
  titles <- c("Genome coverage: expected vs observed", "Target read pairs: designed vs counted")
  xlabels <- c("Expected depth from target FASTQ (x)", "Designed target read pairs")
  ylabels <- c("Observed genome mean depth (x)", "Counted target FASTQ pairs")
  for (i in 1:2) {
    x <- d[[fields[[i]][1]]]
    y <- d[[fields[[i]][2]]]
    valid <- is.finite(x) & is.finite(y)
    upper <- max(c(1, x[valid], y[valid])) * 1.12
    plot(NA, xlim = c(0, upper), ylim = c(0, upper), bty = "l", xlab = xlabels[i], ylab = ylabels[i],
         main = titles[i], cex.main = 0.92)
    abline(a = 0, b = 1, lty = 2, col = "#777777")
    points(x[valid], y[valid], pch = 16, cex = 0.9, col = adjustcolor("#287C8E", alpha.f = 0.65))
    mtext(sprintf("%d evaluable / %d planned runs", sum(valid), nrow(d)), side = 3, line = 0.2, cex = 0.8)
  }
  mtext("RRMG | target-only read validation", side = 3, outer = TRUE, line = 1.1, font = 2, cex = 1.1)
  mtext("Dashed line: expected equality. Actual FASTQ-derived coverage is distinct from designed target pairs.",
        side = 1, outer = TRUE, line = 0.9, cex = 0.8)
  mtext("Overlapping points are retained; inspect summary_source_data.tsv for per-run counts and ratios.",
        side = 1, outer = TRUE, line = 2, cex = 0.8)
}

save_design_figures <- function(bundle, out) {
  d <- bundle$data
  summary <- design_summary(d)
  write_tsv(summary, file.path(out, "detection_by_design.tsv"))
  save_pair(file.path(out, "design_hit_rates"), 11, max(5, 3 + 0.65 * length(unique(d$TargetAbundanceDesign))),
            function() draw_hit_rates(summary))
  sets <- list(
    contig_metrics_by_design = list(title = "Target-gene assembly recovery by design",
      metrics = c(ContigGeneBreadthPct = "Gene breadth (%)", ContigCoveredBp = "Covered bases (bp)",
                  ContigBestAlnLen = "Best alignment (bp)", ContigBestIdentityPct = "Best identity (%)")),
    read_metrics_by_design = list(title = "Mixed QC reads: target-gene mapping by design",
      metrics = c(ReadMappedAlignments = "Mapped alignments", ReadGeneBreadthPct = "Gene breadth (%)",
                  ReadGeneMeanDepth = "Mean depth (x)", ReadMaxDepth = "Max depth (x)", ReadBreadthGE5Pct = "Breadth at depth >=5 (%)")),
    target_genome_metrics_by_design = list(title = "Target-only genome coverage by design",
      metrics = c(TargetCovFromFastq = "FASTQ expected depth (x)", TargetGenomeMeanDepth = "Observed mean depth (x)",
                  TargetGenomeBreadthPct = "Genome breadth (%)", GenomeDepthObservedOverExpected = "Observed/expected depth"))
  )
  abundance <- sort(unique(d$TargetAbundanceDesign))
  pages <- split(abundance, ceiling(seq_along(abundance) / 3))
  for (name in names(sets)) {
    for (p in seq_along(pages)) {
      save_pair(file.path(out, sprintf("%s_%02d", name, p)), max(7, 4 * length(pages[[p]])),
                1.5 + 2.3 * length(sets[[name]]$metrics),
                function() draw_grouped_metrics(d, sets[[name]]$metrics, pages[[p]], sets[[name]]$title))
    }
  }
  save_pair(file.path(out, "target_expected_vs_observed"), 11, 5.4, function() draw_validation(d))
}

save_pair <- function(stem, width, height, draw) {
  for (format in c("png", "pdf")) {
    path <- paste0(stem, ".", format)
    if (format == "png") {
      png(path, width = width, height = height, units = "in", res = 180, type = "cairo", bg = "white")
    } else {
      cairo_pdf(path, width = width, height = height, family = "sans", bg = "white")
    }
    tryCatch(draw(), finally = dev.off())
  }
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2 || any(args %in% c("-h", "--help"))) {
    cat("Usage: Rscript rrmg_plot_preview.R WORKDIR OUTPUT_DIR\n",
        "Reads collected tables + gene depths; writes PNG/PDF previews and source/audit TSVs.\n",
        "Managed outputs are replaced after successful rendering; no analysis input is changed.\n", sep = "")
    if (length(args) == 1 && args[1] %in% c("-h", "--help")) return(invisible(NULL))
    stop("Expected exactly two paths.")
  }
  wd <- normalizePath(args[1], mustWork = TRUE)
  out <- args[2]
  if (!isTRUE(capabilities("cairo")) || !isTRUE(capabilities("png"))) {
    stop("R needs cairo + PNG support for headless export; no outputs written.")
  }
  bundle <- collect_preview(wd)
  job <- rrmg_output_begin(out, "rrmg_summary_figures")
  on.exit(rrmg_output_abort(job), add = TRUE)
  out <- job$stage
  write_tsv(bundle$data, file.path(out, "summary_source_data.tsv"))
  write_tsv(bundle$status, file.path(out, "sample_status.tsv"))
  write_tsv(bundle$depth, file.path(out, "gene_depth_source_data.tsv"))
  checksums <- tools::md5sum(bundle$inputs)
  write_tsv(data.frame(File = names(checksums), MD5 = unname(checksums)), file.path(out, "input_checksums.tsv"))
  writeLines(if (length(bundle$notes)) bundle$notes else "No preview checks raised warnings.", file.path(out, "warnings.txt"))
  capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
  save_design_figures(bundle, out)
  overview_pages <- split(seq_len(nrow(bundle$data)), ceiling(seq_len(nrow(bundle$data)) / 12))
  for (p in seq_along(overview_pages)) {
    save_pair(file.path(out, sprintf("overview_%02d", p)), 10, 8,
              function() draw_overview(bundle, overview_pages[[p]], p, length(overview_pages)))
  }
  depth_pages <- split(seq_along(bundle$panels), ceiling(seq_along(bundle$panels) / 3))
  for (p in seq_along(depth_pages)) {
    panels <- bundle$panels[depth_pages[[p]]]
    save_pair(file.path(out, sprintf("gene_depth_%02d", p)), 10, 1.5 + 2.6 * length(panels),
              function() draw_depth(bundle, panels, p, length(depth_pages)))
  }
  writeLines(c(
    "RRMG local output previews",
    "Adapted figure families: RRMG_permutation_summary_plots.R and RRMG_Afum_target_reads_genome_coverage_plots.R.",
    "Implementation: portable base R; the original ggplot2 appearance is not reproduced exactly.",
    "design_hit_rates: descriptive contig/read hit fractions. Denominator = evaluable runs, with planned n shown.",
    "detection_by_design.tsv: hits, valid denominators and missing counts for contig/read/both statuses.",
    "*_metrics_by_design: original metric families, grouped by designed abundance and total pairs; boxes require n>=3.",
    "target_expected_vs_observed: genome depth from actual target FASTQ versus mapped depth; designed versus actual pair counts.",
    "overview_*.png/pdf: a = contig/mixed-QC-read gene breadth; b = mixed-QC-read gene mean depth;",
    "c = target-only genome breadth; d = designed versus counted target-only FASTQ pairs.",
    "gene_depth_*.png/pdf: raw per-base depth, one sample/reference record per panel, no smoothing.",
    "Each point is one simulation run, not a biological replicate. Descriptive groups only; no statistical tests.",
    "Row labels map to RunID in summary_source_data.tsv; plan order is retained.",
    "Summary values come directly from collected tables; no NA-to-zero conversion or imputation.",
    "Sample status reports presence and selected consistency checks, not full biological validation.",
    "Both summaries and depth files must describe the same analysis run; inspect warnings.txt.",
    "Target-only genome coverage is not mixed-metagenome species abundance. No Kraken/VAMB/AMR conclusions.",
    "VarN is not plotted: the upstream placeholder is not a verified mutation call.",
    "Alignment-block and FungAMR-site overlays: use the subsequent input preparation + alignment plotting sections of runbook.sh.",
    "Target genome window plots are not generated by this summary helper.",
    "PNG: 180 dpi preview. PDF: vector plots with selectable text. Input checksums and R versions included."
  ), file.path(out, "README.txt"))
  rrmg_output_finish(job)
  cat(sprintf("[FIGURES] %s\n[RUNS] %d planned; %d preview warning(s)\n", job$final, nrow(bundle$data), length(bundle$notes)))
  if (length(bundle$notes)) cat("[CHECK] Read warnings.txt before interpreting the plots.\n")
  cat("[DONE] PNG + PDF + source data; input analysis files were not modified.\n")
  invisible(bundle)
}

if (sys.nframe() == 0L) main()
