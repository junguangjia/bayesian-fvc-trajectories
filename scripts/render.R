#!/usr/bin/env Rscript
source("scripts/common.R")
args <- commandArgs(trailingOnly = TRUE)
run_id <- arg_value(args, "--run", "demo")
if (!grepl("^[a-zA-Z0-9_-]+$", run_id)) stop("Invalid run identifier.")
run_dir <- normalizePath(file.path("results", run_id), mustWork = TRUE)
manifest <- jsonlite::read_json(file.path(run_dir, "manifest.json"), simplifyVector = TRUE)
if (!manifest$status %in% c("verified", "demo") ||
    is.null(manifest$completed_utc) || manifest$run_count != manifest$planned_run_count ||
    (manifest$mode == "research" && manifest$status != "verified") ||
    (manifest$mode != "research" && manifest$status != "demo")) stop("Cannot render incomplete or unverified results.")
if (is.null(manifest$output_sha256)) stop("Run artifact hashes are missing.")
for (f in names(manifest$output_sha256)) {
  target <- file.path(run_dir, f)
  if (!file.exists(target) || digest::digest(file = target, algo = "sha256") != manifest$output_sha256[[f]])
    stop("Run artifact integrity check failed.")
}
source("scripts/render_latex.R")
if (!"--tex-only" %in% args) {
  engine <- Sys.which("pdflatex")
  if (!nzchar(engine)) stop("LaTeX source created; pdflatex unavailable. Install a LaTeX distribution or use --tex-only and the built-in LaTeX editor.")
  build <- file.path("private", "latex", run_id)
  dir.create(build, recursive = TRUE, showWarnings = FALSE)
  for (pass in 1:2) {
    status <- system2(engine, c("-interaction=nonstopmode", "-halt-on-error",
      shQuote(paste0("-output-directory=", build)), shQuote(file.path(run_dir, "report.tex"))),
      stdout = file.path(build, "compile.txt"), stderr = file.path(build, "compile-errors.txt"))
    if (status != 0L) stop("LaTeX compilation failed; inspect the private build log.")
  }
  if (!file.copy(file.path(build, "report.pdf"), file.path(run_dir, "report.pdf"), overwrite = TRUE)) stop("Could not export report PDF.")
  cat("PDF:", file.path("results", run_id, "report.pdf"), "\n")
}
