#!/usr/bin/env Rscript
# Create an aggregate-only validation report; do not compile or publish it.
# Receipt schema (strict booleans; both hashes must match the actual manifests):
# {"schema_version":1,"research_analysis_sha256":"...",
#  "recovery_analysis_sha256":"...","unit_numerics_passed":true,
#  "compiled_smoke_passed":true,"clean_reproduction_passed":true}
args <- commandArgs(trailingOnly = TRUE)
evidence_index <- match("--evidence", args)
evidence_path <- if (!is.na(evidence_index) && evidence_index < length(args)) args[evidence_index + 1L] else NULL
research_dir <- file.path("results", "research")
recovery_dir <- file.path("results", "recovery")
if (!dir.exists(research_dir)) stop("The research run directory does not exist.")
read_json <- function(path, simplify = TRUE) {
  if (!file.exists(path)) return(list())
  tryCatch(jsonlite::fromJSON(path, simplifyVector = simplify), error = function(e) list())
}
read_public <- function(path, required_columns) {
  if (!file.exists(path)) return(data.frame())
  out <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) data.frame())
  if (!all(required_columns %in% names(out))) data.frame() else out
}
research <- read_json(file.path(research_dir, "manifest.json"))
recovery <- read_json(file.path(recovery_dir, "manifest.json"))
receipt <- if (is.null(evidence_path)) list() else read_json(evidence_path)
diagnostic_columns <- c("run", "model", "variant", "max_rhat", "min_bulk_ess", "min_tail_ess",
  "max_slope_mcse", "divergences", "treedepth_hits", "min_ebfmi", "passed")
research_diagnostics <- read_public(file.path(research_dir, "diagnostics.csv"), diagnostic_columns)
recovery_diagnostics <- read_public(file.path(recovery_dir, "diagnostics.csv"), diagnostic_columns)
recovery_summary <- read_public(file.path(recovery_dir, "recovery.csv"),
  c("estimand", "truth", "estimate", "absolute_error", "tolerance_ml_week", "passed"))
checks <- read_public(file.path(research_dir, "checks.csv"),
  c("model", "statistic", "observed", "lower", "upper", "ppp"))
metadata <- read_json(file.path(research_dir, "fit_metadata.json"), simplify = FALSE)

scalar <- function(x, fallback = "not recorded") {
  if (is.null(x) || length(x) != 1L || is.na(x)) fallback else as.character(x)
}
fmt <- function(x, digits = 3L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x)) return("not available")
  formatC(x, format = "f", digits = digits)
}
tex <- function(x) {
  escapes <- c("\\" = "\\textbackslash{}", "{" = "\\{", "}" = "\\}", "#" = "\\#",
    "$" = "\\$", "%" = "\\%", "&" = "\\&", "_" = "\\_", "^" = "\\textasciicircum{}", "~" = "\\textasciitilde{}")
  vapply(as.character(x), function(value) {
    if (is.na(value)) return("not available")
    chars <- strsplit(enc2utf8(value), "", fixed = TRUE)[[1L]]
    replacement <- unname(escapes[chars])
    replacement[is.na(replacement)] <- chars[is.na(replacement)]
    paste0(replacement, collapse = "")
  }, character(1L), USE.NAMES = FALSE)
}
true <- function(x) !is.na(x) & tolower(as.character(x)) %in% c("true", "1")
finite_column <- function(x) is.numeric(x) && length(x) > 0L && all(is.finite(x))
same_hash <- function(a, b) is.character(a) && length(a) == 1L && !is.na(a) &&
  grepl("^[a-f0-9]{64}$", a) && identical(a, b)
manifest_ok <- function(manifest, mode, status, count) {
  identical(manifest$mode, mode) && identical(manifest$status, status) &&
    is.numeric(manifest$run_count) && length(manifest$run_count) == 1L &&
    isTRUE(manifest$run_count == count) && is.numeric(manifest$seed) &&
    length(manifest$seed) == 1L && isTRUE(manifest$seed == 4224)
}
diagnostics_ok <- function(x, expected) {
  if (!nrow(x) || nrow(x) != length(expected)) return(FALSE)
  keys <- paste(x$model, x$variant, sep = "_")
  if (anyDuplicated(keys) || !setequal(keys, expected) || !all(true(x$passed))) return(FALSE)
  numeric_columns <- c("max_rhat", "min_bulk_ess", "min_tail_ess", "max_slope_mcse", "divergences", "treedepth_hits", "min_ebfmi")
  if (!all(vapply(x[numeric_columns], finite_column, logical(1L)))) return(FALSE)
  all(x$max_rhat < 1.01 & x$min_bulk_ess >= 400 & x$min_tail_ess >= 400 &
        x$divergences == 0 & x$treedepth_hits == 0 & x$min_ebfmi >= .3 & x$max_slope_mcse >= 0)
}
artifact_hashes_ok <- function(manifest, directory, files) {
  if (is.null(manifest$output_sha256)) return(FALSE)
  all(vapply(files, function(name) {
    path <- file.path(directory, name)
    expected <- manifest$output_sha256[[name]]
    file.exists(path) && same_hash(expected, digest::digest(file = path, algo = "sha256"))
  }, logical(1L)))
}
research_expected <- c(paste0(c("P", "H", "C"), "_primary"),
  as.vector(outer(c("P", "H", "C"), paste0("_cv_fold_", 1:5), paste0)),
  as.vector(outer(c("H", "C"), c("_prior_half", "_prior_double", "_all_rows"), paste0)))
recovery_expected <- paste0(c("P", "H", "C"), "_primary")
meta_keys <- vapply(metadata, function(entry) paste(scalar(entry$model), scalar(entry$variant), sep = "_"), character(1L))
attempt_counts <- vapply(metadata, function(entry) length(entry$attempts), integer(1L))
attempts <- if (length(attempt_counts)) sum(attempt_counts) else NA_integer_
attempts_ok <- length(metadata) == 24L && !anyDuplicated(meta_keys) &&
  setequal(meta_keys, research_expected) && all(attempt_counts %in% 1:3) &&
  is.numeric(research$sampling_attempts) && length(research$sampling_attempts) == 1L &&
  isTRUE(research$sampling_attempts == attempts)
research_hashes <- artifact_hashes_ok(research, research_dir,
  c("diagnostics.csv", "fit_metadata.json", "checks.csv"))
recovery_hashes <- artifact_hashes_ok(recovery, recovery_dir, c("diagnostics.csv", "recovery.csv"))
recovery_values_ok <- nrow(recovery_summary) == 1L &&
  identical(recovery_summary$estimand, "standardized_slope") &&
  all(vapply(recovery_summary[c("truth", "estimate", "absolute_error", "tolerance_ml_week")], finite_column, logical(1L))) &&
  isTRUE(recovery_summary$tolerance_ml_week > 0) &&
  isTRUE(abs(recovery_summary$absolute_error - abs(recovery_summary$estimate - recovery_summary$truth)) < 1e-8) &&
  isTRUE(recovery_summary$absolute_error < recovery_summary$tolerance_ml_week) && all(true(recovery_summary$passed))
receipt_schema_ok <- is.numeric(receipt$schema_version) && length(receipt$schema_version) == 1L && isTRUE(receipt$schema_version == 1)
receipt_bound <- receipt_schema_ok && same_hash(receipt$research_analysis_sha256, research$analysis_sha256) &&
  same_hash(receipt$recovery_analysis_sha256, recovery$analysis_sha256)
gates <- c(
  "Verified research manifest and seed" = manifest_ok(research, "research", "verified", 24L),
  "All 24 research configurations and diagnostic gates" = diagnostics_ok(research_diagnostics, research_expected),
  "Attempt counts agree with the manifest" = attempts_ok,
  "Research artifact hashes match" = research_hashes,
  "Synthetic recovery manifest and seed" = manifest_ok(recovery, "recovery", "demo", 3L),
  "All three recovery models meet diagnostic gates" = diagnostics_ok(recovery_diagnostics, recovery_expected),
  "Recovery error is below its recorded tolerance" = recovery_values_ok,
  "Recovery artifact hashes match" = recovery_hashes,
  "QA receipt matches both analysis hashes" = receipt_bound,
  "Unit and numerical tests" = receipt_bound && isTRUE(receipt$unit_numerics_passed),
  "Compiled synthetic smoke workflow" = receipt_bound && isTRUE(receipt$compiled_smoke_passed),
  "Clean-directory synthetic reproduction" = receipt_bound && isTRUE(receipt$clean_reproduction_passed))
complete <- all(gates)
serial <- if (nrow(checks)) checks[checks$model %in% c("H", "C") & checks$statistic == "residual_lag1_correlation", , drop = FALSE] else data.frame()
serial_flags <- if (nrow(serial)) serial[is.finite(serial$ppp) & (serial$ppp < .025 | serial$ppp > .975), , drop = FALSE] else data.frame()

lines <- character()
emit <- function(...) lines <<- c(lines, paste0(...))
paragraph <- function(...) { emit(...); emit("") }
row <- function(cells) emit(paste(cells, collapse = " & "), " \\\\")
table_start <- function(spec, headers) {
  spec <- gsub("p{", ">{\\raggedright\\arraybackslash}p{", spec, fixed = TRUE)
  emit("\\begingroup\\small\\setlength{\\tabcolsep}{5pt}\\begin{longtable}{", spec, "}\\toprule")
  row(tex(headers)); emit("\\midrule\\endfirsthead\\toprule")
  row(tex(headers)); emit("\\midrule\\endhead\\bottomrule\\endfoot")
}
table_end <- function() emit("\\end{longtable}\\endgroup", "")
extreme <- function(x, column, fn) {
  if (!nrow(x) || !finite_column(x[[column]])) NA_real_ else fn(x[[column]])
}
emit(c("\\documentclass[11pt,letterpaper]{article}", "\\usepackage[margin=1in]{geometry}",
  "\\usepackage[T1]{fontenc}", "\\usepackage{lmodern}", "\\usepackage{amsmath,booktabs,longtable,array}",
  "\\usepackage[hidelinks]{hyperref}", "\\hypersetup{pdftitle={Computational Validation of Bayesian FVC Trajectory Models},pdfauthor={Junguang Jia}}", "\\setlength{\\parindent}{0pt}", "\\setlength{\\parskip}{0.3em}",
  "\\setlength{\\emergencystretch}{2em}", "\\title{Computational Validation of Bayesian FVC Trajectory Models}",
  "\\author{Junguang Jia}", "\\date{}", "\\begin{document}", "\\maketitle"))
emit("\\section{Scope and disposition}")
paragraph("This report evaluates implementation checks, sampling diagnostics, and a single synthetic parameter-recovery experiment. It reads the recorded research and recovery outputs and an explicit QA receipt. The common base seed is 4224. It does not certify external predictive validity, clinical utility, causal interpretation, or general parameter-recovery calibration.")
paragraph("\\textbf{Disposition: ", if (complete) "computational verification complete." else "incomplete; verification is not established.", "} ",
  if (nrow(serial_flags)) "Residual checks additionally identify a substantive serial-dependence limitation, independently of computational acceptance." else
    "Model adequacy remains a separate question from numerical verification.")

emit("\\section{Research configuration and sampling checks}")
paragraph("The prespecified research analysis has 24 configurations: three main fits, 15 patient-fold fits, four prior-scale variants, and two duplicate-handling variants. There are ",
  fmt(nrow(research_diagnostics), 0L), " recorded diagnostic rows and ", fmt(attempts, 0L),
  " recorded sampling attempts", if (attempts_ok) paste0(" (", fmt(attempts - 24L, 0L), " retries)") else "",
  ". Attempts are counted separately from configurations; an initial failed attempt is retained in the metadata.")
table_start("lrrrrrrr", c("Analysis", "Fits", "Max R-hat", "Min bulk", "Min tail", "Div.", "Depth", "Min BFMI"))
for (pair in list(list("Research", research_diagnostics), list("Recovery", recovery_diagnostics))) {
  d <- pair[[2L]]
  row(c(tex(pair[[1L]]), fmt(nrow(d), 0L), fmt(extreme(d, "max_rhat", max), 4L),
    fmt(extreme(d, "min_bulk_ess", min), 0L), fmt(extreme(d, "min_tail_ess", min), 0L),
    fmt(extreme(d, "divergences", sum), 0L), fmt(extreme(d, "treedepth_hits", sum), 0L),
    fmt(extreme(d, "min_ebfmi", min), 3L)))
}
table_end()
paragraph("Gates are rank-normalized R-hat $<1.01$, bulk and tail ESS $\\geq400$, no post-warmup divergences, no tree-depth hits, and per-chain E-BFMI $\\geq0.3$. Research slope MCSE has recorded maximum ",
  fmt(extreme(research_diagnostics, "max_slope_mcse", max), 4L),
  " mL/week. The report recomputes gates from the extrema rather than relying only on the saved pass flag. These numerical gates do not establish that the Gaussian residual model is adequate.")

emit("\\section{Synthetic parameter recovery}")
if (nrow(recovery_summary)) {
  table_start("lrrrr", c("Estimand", "Truth", "Estimate", "Absolute error", "Tolerance"))
  for (i in seq_len(nrow(recovery_summary))) {
    z <- recovery_summary[i, ]
    row(c(if (z$estimand == "standardized_slope") "C overall slope" else tex(z$estimand),
      fmt(z$truth, 4L), fmt(z$estimate, 4L), fmt(z$absolute_error, 4L), fmt(z$tolerance_ml_week, 4L)))
  }
  table_end()
} else paragraph("No valid recovery summary is available.")
paragraph("All quantities in this table are in mL/week. The acceptance rule is absolute error strictly below the recorded tolerance. The calculation is independently checked against the saved truth and estimate. This is one seeded synthetic dataset and one targeted slope-recovery check; it is not simulation-based calibration or a repeated-simulation coverage study.")

emit("\\section{Model-adequacy limitation}")
if (nrow(serial)) {
  table_start("lrrrl", c("Model", "Observed", "Replicated 95% interval", "Tail fraction", "Flag"))
  for (i in seq_len(nrow(serial))) {
    z <- serial[i, ]
    flag <- is.finite(z$ppp) && (z$ppp < .025 || z$ppp > .975)
    row(c(tex(z$model), fmt(z$observed, 3L), paste0("[", fmt(z$lower, 3L), ", ", fmt(z$upper, 3L), "]"),
      fmt(z$ppp, 3L), if (flag) "yes" else "no"))
  }
  table_end()
} else paragraph("The residual serial-dependence check is unavailable.")
paragraph("The statistic is pooled adjacent-visit residual correlation within patients ordered by time. Observed entries are medians of draw-specific discrepancies; tail fractions compare replicated and observed discrepancies within each posterior draw. Fractions below 0.025 or above 0.975 flag lack of agreement, not calibrated hypothesis-test significance. ",
  if (nrow(serial_flags)) "The flagged discrepancy indicates remaining serial dependence. Nominal population credible intervals and new-patient predictive intervals retain a questionable residual-independence assumption, so calibration and predictive adequacy may be affected." else
    "An unflagged check would not itself establish residual independence.")

emit("\\section{Evidence receipt and reproducibility}")
table_start("p{0.76\\textwidth}l", c("Required evidence", "Result"))
for (i in seq_along(gates)) row(c(tex(names(gates)[i]), if (gates[i]) "passed" else "not established"))
table_end()
paragraph("Software checks are asserted only when the supplied receipt contains literal true values and its analysis hashes match both run manifests. Missing, false, malformed, or mismatched evidence cannot establish a pass. Research diagnostic, fit-metadata, and posterior-check hashes and recovery diagnostic and recovery-summary hashes are compared with the recorded artifact hashes. The receipt's local path and host details are not published.")
paragraph("The receipt covers unit and numerical checks, a compiled synthetic smoke workflow, and reproduction from a clean directory. These checks complement the full research fits; a smoke run alone does not establish convergence. This report makes no hosted-CI, deployment, publication, independent-data, or clinical-validation claim.")
emit("\\end{document}")
writeLines(lines, file.path(research_dir, "validation.tex"), useBytes = TRUE)
cat("Validation LaTeX written to results/research/validation.tex. Status:",
  if (complete) "computational verification complete" else "incomplete", "\n")
if (!complete) quit(save = "no", status = 1L)
