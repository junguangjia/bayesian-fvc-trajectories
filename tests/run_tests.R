#!/usr/bin/env Rscript
source("scripts/common.R")
fail_expected <- function(x) stopifnot(inherits(tryCatch({force(x); NULL}, error = identity), "error"))
raw <- simulate_data(15)
d <- prepare_data(raw)
stopifnot(nrow(d$patients) == 15L, nrow(d$obs) == 120L, all(d$obs$y > 0),
  identical(colnames(design_matrix(d$obs, "C")), c("alpha", "beta_time", "beta_age", "beta_male",
    "beta_current", "beta_never", "gamma_age", "gamma_male", "gamma_current", "gamma_never")))
shuffled <- prepare_data(raw[sample(nrow(raw)), ])
stopifnot(identical(d, shuffled), all(d$obs$time == d$obs$week / 52))
bad <- raw; bad$Sex[1] <- "unrecognized"; fail_expected(prepare_data(bad))
bad <- raw; bad$Age[1] <- NA; fail_expected(prepare_data(bad))
bad <- raw; bad$FVC[1] <- Inf; fail_expected(prepare_data(bad))
bad <- raw; bad$Age[1] <- bad$Age[1] + 1; fail_expected(prepare_data(bad))
fail_expected(prepare_data(raw[, -1])); fail_expected(prepare_data(rbind(raw, raw[1, ])))
repeat_row <- raw[1, ]; repeat_row$FVC <- repeat_row$FVC + 100
replicated <- prepare_data(rbind(raw, repeat_row))
stopifnot(replicated$audit$repeated_week_extra_rows == 1,
  abs(replicated$obs$y[1] - d$obs$y[1] - .05) < 1e-12,
  nrow(prepare_data(rbind(raw, repeat_row), FALSE)$obs) == nrow(raw) + 1)
folds <- patient_folds(d)
stopifnot(identical(folds, patient_folds(d)), !anyDuplicated(folds$patient))
for (k in 1:5) {
  train <- subset_prepared(d, folds$patient[folds$fold != k])
  test <- subset_prepared(d, folds$patient[folds$fold == k])
  stopifnot(!length(intersect(train$obs$patient, test$obs$patient)),
    nrow(train$obs) + nrow(test$obs) == nrow(d$obs),
    identical(sort(unique(train$obs$id)), seq_len(nrow(train$patients))))
}
b <- matrix(0, 2, 10); b[, 2] <- -.25; b[, 7:10] <- .1
draws <- list(beta = b)
expected <- (-.25 + .1 * sum(colMeans(d$patients[, c("age10", "male", "current", "never")]))) * 1000 / 52
stopifnot(all(abs(slope_draws(draws, d, "C") - expected) < 1e-12))
source("tests/evaluate_tests.R")
run_evaluate_tests()
# Run the model's input/diagnostic tests in an isolated environment: its stub
# design matrix must not mask the real preprocessing function above.
source("tests/test-model.R", local = new.env(parent = globalenv()))
cat("PASS: data, design, units, duplicate handling, grouped folds, derived slopes, prediction and diagnostics.\n")
