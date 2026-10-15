# Managed figure directories: render first, publish only after success.
# The inventory prevents replacement of unrelated or manually edited files.
rrmg_output_inventory <- function(path) {
  files <- list.files(path, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  files <- setdiff(files, ".rrmg_outputs.tsv")
  links <- list.files(path, recursive = TRUE, all.files = TRUE, no.. = TRUE, include.dirs = TRUE)
  if (any(nzchar(Sys.readlink(c(path, file.path(path, links)))), na.rm = TRUE))
    stop("Output symlinks are not supported: ", path)
  data.frame(File = files, MD5 = unname(tools::md5sum(file.path(path, files))), stringsAsFactors = FALSE)
}

rrmg_output_check <- function(path, kind) {
  if (!file.exists(path)) return(invisible(NULL))
  marker <- file.path(path, ".rrmg_outputs.tsv")
  if (!dir.exists(path) || !file.exists(marker))
    stop("Existing output is not a managed ", kind, " directory; move it aside first: ", path)
  old <- read.delim(marker, stringsAsFactors = FALSE, check.names = FALSE)
  now <- rrmg_output_inventory(path)
  if (!all(c("Kind", "File", "MD5") %in% names(old)) || !nrow(old) ||
      anyNA(old) || any(old$Kind != kind) || !identical(old$File, now$File) ||
      !identical(old$MD5, now$MD5))
    stop("Output contains added, changed or missing files; preserve/move it aside first: ", path)
}

rrmg_output_begin <- function(path, kind) {
  parent <- dirname(path)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
  final <- file.path(parent, basename(path))
  lock <- paste0(final, ".lock")
  if (!dir.create(lock, showWarnings = FALSE))
    stop("Output is locked: ", lock, "; check for another running process before removing the empty lock directory.")
  ok <- FALSE
  on.exit(if (!ok) unlink(lock, recursive = TRUE), add = TRUE)
  rrmg_output_check(final, kind)
  stage <- tempfile(paste0(".", basename(final), "-stage-"), tmpdir = parent)
  if (!dir.create(stage)) stop("Cannot create staging directory: ", stage)
  ok <- TRUE
  list(final = final, stage = stage, lock = lock, kind = kind)
}

rrmg_output_abort <- function(job) {
  # Both paths were created by rrmg_output_begin, never supplied input directories.
  if (dir.exists(job$stage)) unlink(job$stage, recursive = TRUE)
  if (dir.exists(job$lock)) unlink(job$lock, recursive = TRUE)
}

rrmg_output_finish <- function(job) {
  inventory <- rrmg_output_inventory(job$stage)
  if (!nrow(inventory)) stop("No figure outputs generated")
  inventory <- cbind(Kind = job$kind, inventory)
  write.table(inventory, file.path(job$stage, ".rrmg_outputs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  rrmg_output_check(job$final, job$kind)
  backup <- tempfile(paste0(".", basename(job$final), "-previous-"), tmpdir = dirname(job$final))
  had_old <- dir.exists(job$final)
  if (had_old && !file.rename(job$final, backup)) stop("Cannot stage previous outputs: ", job$final)
  if (!file.rename(job$stage, job$final)) {
    if (had_old && !file.rename(backup, job$final)) warning("Previous outputs remain at: ", backup)
    stop("Cannot publish outputs: ", job$final)
  }
  # Only remove the verified, inventoried previous directory after publication.
  if (had_old) unlink(backup, recursive = TRUE)
  rrmg_output_abort(job)
  invisible(job$final)
}
