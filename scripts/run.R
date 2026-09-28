#!/usr/bin/env Rscript
source("scripts/common.R")
args <- commandArgs(trailingOnly = TRUE)
mode <- arg_value(args, "--mode", "demo")
if (!mode %in% c("demo", "research", "recovery")) stop("Unknown mode.")
smoke <- "--smoke" %in% args
if (smoke && mode == "research") stop("Smoke settings cannot generate research results.")
run_id <- arg_value(args, "--run", mode)
if (!grepl("^[a-zA-Z0-9_-]+$", run_id)) stop("Invalid run identifier.")
dir.create(file.path("private", run_id), recursive = TRUE, showWarnings = FALSE, mode = "0700")
dir.create(file.path("results", run_id), recursive = TRUE, showWarnings = FALSE)
dir.create("cache", showWarnings = FALSE, mode = "0700")
Sys.chmod(c("private", file.path("private", run_id), "cache"), mode = "0700")
out_dir <- file.path("results", run_id)
private_dir <- file.path("private", run_id)
# Invalidate an earlier release before touching any of this run's artifacts.
jsonlite::write_json(list(run_id = run_id, mode = mode, status = "unverified", run_count = 0),
  file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
for (artifact in c("report.html", "report.tex", "report.pdf", "validation.tex", "validation.pdf")) {
  target <- file.path(out_dir, artifact)
  if (file.exists(target)) unlink(target)
}
seed <- 4224L
if (mode == "research") {
  path <- arg_value(args, "--data")
  if (is.null(path) || !file.exists(path)) stop("Research mode requires --data <csv>.")
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  data_hash <- digest::digest(file = path, algo = "sha256")
} else {
  raw <- simulate_data(if (smoke) 15L else if (mode == "recovery") 100L else 50L, seed)
  data_hash <- digest::digest(raw, algo = "sha256")
}
prepared <- prepare_data(raw)
original <- prepare_data(raw, aggregate = FALSE)
settings <- if (smoke) list(chains = 2L, iter = 200L, warmup = 100L,
  adapt_delta = .95, max_treedepth = 12L, cores = 2L) else
  list(chains = 4L, iter = 6000L, warmup = 3000L,
       adapt_delta = .95, max_treedepth = 12L, cores = 4L)
write_csv(data.frame(metric = names(prepared$audit), value = unlist(prepared$audit)), file.path(out_dir, "cohort.csv"))
cat("Valid input:", prepared$audit$patients, "patients;", prepared$audit$analyzed_rows, "analysis records.\n")
prior <- do.call(rbind, lapply(c("P", "H", "C"), function(m) prior_predictive(prepared, m, seed)))
write_csv(prior, file.path(out_dir, "prior_checks.csv"))
cat("Prior predictive checks completed before fitting.\n")
jobs <- lapply(c("P", "H", "C"), function(m) list(model = m, variant = "primary", scale = 1, data = prepared))
if (mode == "research") {
  folds <- patient_folds(prepared, seed = seed)
  write_csv(folds, file.path(private_dir, "folds.csv"))
  for (fold in 1:5) for (m in c("P", "H", "C")) {
    train <- subset_prepared(prepared, folds$patient[folds$fold != fold])
    test <- subset_prepared(prepared, folds$patient[folds$fold == fold])
    jobs[[length(jobs) + 1L]] <- list(model = m, variant = paste0("cv_fold_", fold), scale = 1,
                                     data = train, test = test, fold = fold)
  }
  for (scale in c(.5, 2)) for (m in c("H", "C"))
    jobs[[length(jobs) + 1L]] <- list(model = m, variant = if (scale == .5) "prior_half" else "prior_double",
                                     scale = scale, data = prepared)
  for (m in c("H", "C")) jobs[[length(jobs) + 1L]] <- list(model = m, variant = "all_rows", scale = 1, data = original)
}
estimates <- diags <- checks <- curves <- cv <- NULL
fit_meta <- list()
software <- as.list(vapply(c("rstan", "StanHeaders", "posterior", "digest", "rmarkdown", "knitr", "jsonlite"),
  function(p) as.character(packageVersion(p)), character(1)))
software$R <- as.character(getRversion())
source_files <- c(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)),
  sort(list.files("stan", pattern = "\\.stan$", full.names = TRUE)), "scripts/run.R", "scripts/common.R", "renv.lock")
source_hashes <- as.list(setNames(vapply(source_files, function(f) digest::digest(file = f, algo = "sha256"), character(1)), source_files))
manifest <- list(mode = if (smoke) "smoke" else mode, run_id = run_id,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), seed = seed, status = "unverified",
  run_count = 0L, planned_run_count = length(jobs), data_sha256 = data_hash, software = software,
  settings = settings, prediction_draw_limit = 1000L,
  source_sha256 = source_hashes, analysis_sha256 = digest::digest(source_hashes, algo = "sha256"))
write_manifest <- function() jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
write_manifest()

for (i in seq_along(jobs)) {
  job <- jobs[[i]]; label <- paste(job$model, job$variant, sep = "_")
  cat(sprintf("[%d/%d] %s\n", i, length(jobs), label)); flush.console()
  actual <- settings
  attempts <- list()
  for (attempt in 1:3) {
    res <- fit_model(job$data, job$model, job$scale, settings = actual, cache_dir = "cache", seed = seed + i)
    d <- diagnostic_summary(res$fit, job$data, job$model, job$variant, label)
    attempts[[attempt]] <- list(settings = actual, diagnostics = d,
      elapsed_seconds = res$metadata$elapsed_seconds, cache_key = res$metadata$content_sha256)
    if (smoke || isTRUE(d$passed)) break
    if (attempt == 1L) {
      cat("Diagnostic gate failed; retaining evidence and retrying with longer chains.\n")
      actual$iter <- 10000L; actual$warmup <- 5000L; actual$adapt_delta <- .99
      if (d$treedepth_hits > 0) actual$max_treedepth <- 14L
    } else if (attempt == 2L) {
      cat("Retry gate remains unmet; applying a targeted final computational extension.\n")
      if (d$treedepth_hits > 0) actual$max_treedepth <- 14L
      if (d$max_rhat >= 1.01 || d$min_bulk_ess < 400 || d$min_tail_ess < 400) {
        actual$iter <- 20000L; actual$warmup <- 5000L
      }
    }
  }
  diags <- rbind(diags, d)
  write_csv(diags, file.path(out_dir, "diagnostics.csv"))
  saveRDS(attempts, file.path(private_dir, paste0(label, "_attempts.rds")))
  fit_meta[[label]] <- list(label = label, model = job$model, variant = job$variant,
    settings = actual, seed = seed + i, metadata = res$metadata, attempts = attempts)
  # Export only a selected, path-free metadata interface from the model layer.
  safe_meta <- lapply(fit_meta, function(m) {
    meta <- m$metadata
    list(label = m$label, model = m$model, variant = m$variant, settings = m$settings,
         seed = m$seed, elapsed_seconds = meta$elapsed_seconds, cache_key = meta$content_sha256,
         source_sha256 = meta$source_sha256, implementation_sha256 = meta$implementation_sha256,
         software_sha256 = meta$software_sha256, settings_sha256 = meta$settings_sha256,
         attempts = m$attempts)
  })
  jsonlite::write_json(safe_meta, file.path(out_dir, "fit_metadata.json"), auto_unbox = TRUE, pretty = TRUE)
  if (!smoke && !isTRUE(d$passed)) {
    manifest$run_count <- i; write_manifest()
    stop("Diagnostic gate remains unmet. Partial results are unverified; inspect private diagnostics.")
  }
  draws <- extract_parameter_draws(res$fit)
  if (grepl("^cv_fold_", job$variant)) {
    ev <- evaluate_new_patients(draws, job$test, job$model, seed = seed + 100L + i)
    cv <- rbind(cv, cbind(model = job$model, fold = job$fold, ev$per_patient))
    write_csv(cv, file.path(private_dir, "cv_patient_scores.csv"))
  } else {
    estimates <- rbind(estimates, estimate_summary(draws, job$data, job$model, job$variant))
    write_csv(estimates, file.path(out_dir, "estimates.csv"))
    if (job$variant == "primary") {
      check <- posterior_checks(draws, job$data, job$model, seed = seed + 200L + i)
      checks <- rbind(checks, cbind(model = job$model, check))
      curves <- rbind(curves, population_curves(draws, job$data, job$model))
      write_csv(checks, file.path(out_dir, "checks.csv")); write_csv(curves, file.path(out_dir, "curves.csv"))
    }
  }
  manifest$run_count <- i; write_manifest()
  rm(draws, res); gc(verbose = FALSE)
}
if (!is.null(cv)) {
  scores <- aggregate_cv(cv)
  write_csv(scores$metrics, file.path(out_dir, "cv_metrics.csv"))
  write_csv(scores$comparisons, file.path(out_dir, "cv_comparisons.csv"))
} else {
  write_csv(data.frame(model = character(), mean_log_score = numeric(), rmse = numeric(), mae = numeric(),
    coverage95 = numeric(), width95 = numeric()), file.path(out_dir, "cv_metrics.csv"))
  write_csv(data.frame(model_a = character(), model_b = character(), difference = numeric(), se = numeric()),
    file.path(out_dir, "cv_comparisons.csv"))
}
if (mode == "recovery") {
  truth <- mean(-.25 + .03 * prepared$patients$age10 - .04 * prepared$patients$male) * 1000 / 52
  c_slope <- estimates[estimates$model == "C" & estimates$estimand == "standardized_slope", ]
  recovery <- data.frame(estimand = "standardized_slope", truth = truth, estimate = c_slope$mean,
    absolute_error = abs(c_slope$mean - truth), tolerance_ml_week = 1.5,
    passed = abs(c_slope$mean - truth) < 1.5)
  write_csv(recovery, file.path(out_dir, "recovery.csv"))
  if (!all(recovery$passed)) stop("Synthetic slope recovery exceeded the prespecified 1.5 mL/week tolerance.")
}
manifest$status <- if (mode == "research") "verified" else "demo"
manifest$completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
manifest$sampling_attempts <- sum(vapply(fit_meta, function(m) length(m$attempts), integer(1)))
artifact_files <- c(list.files(out_dir, pattern = "\\.csv$"), "fit_metadata.json")
manifest$output_sha256 <- as.list(setNames(vapply(artifact_files, function(f)
  digest::digest(file = file.path(out_dir, f), algo = "sha256"), character(1)), artifact_files))
write_manifest()
cat("Completed", length(jobs), "fits. Results:", out_dir, "Status:", manifest$status, "\n")
