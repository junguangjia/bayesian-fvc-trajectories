slope_draws <- function(draws, prepared, model) {
  out <- draws$beta[, 2]
  if (model == "C") out <- out + as.vector(draws$beta[, 7:10, drop = FALSE] %*%
      colMeans(prepared$patients[, c("age10", "male", "current", "never")]))
  out * 1000 / 52
}

summary_row <- function(x, estimand) {
  data.frame(estimand = estimand, mean = mean(x), median = median(x),
    lower95 = unname(quantile(x, .025)), upper95 = unname(quantile(x, .975)),
    p_negative = mean(x < 0))
}

estimate_summary <- function(draws, prepared, model, variant) {
  out <- summary_row(slope_draws(draws, prepared, model), "standardized_slope")
  if (model == "C") {
    out <- rbind(out, summary_row(draws$beta[, 2] * 1000 / 52, "reference_slope"))
    labels <- c("age_main_ml_per_decade", "male_main_ml", "current_main_ml", "never_main_ml",
      "age_time_ml_week_per_decade", "male_time_ml_week", "current_time_ml_week", "never_time_ml_week")
    for (k in 3:10) out <- rbind(out, summary_row(draws$beta[, k] * if (k < 7) 1000 else 1000 / 52, labels[k - 2]))
    profile_columns <- c(male_ex_smoker_age65_slope = 8, female_never_smoker_age65_slope = 10,
      female_current_smoker_age65_slope = 9, female_ex_smoker_age75_slope = 7)
    for (profile in names(profile_columns)) out <- rbind(out,
      summary_row((draws$beta[, 2] + draws$beta[, profile_columns[[profile]]]) * 1000 / 52, profile))
  }
  if (model != "P") {
    out <- rbind(out, summary_row(draws$tau[, 1] * 1000, "tau_intercept"),
      summary_row(draws$tau[, 2] * 1000 / 52, "tau_slope"), summary_row(draws$rho, "rho"))
  }
  cbind(model = model, variant = variant, out)
}

derived_diagnostics <- function(fit, prepared, model) {
  a <- rstan::extract(fit, pars = "beta", permuted = FALSE, inc_warmup = FALSE)
  slope <- a[, , 2, drop = TRUE]
  if (model == "C") {
    means <- colMeans(prepared$patients[, c("age10", "male", "current", "never")])
    for (k in 1:4) slope <- slope + a[, , k + 6, drop = TRUE] * means[k]
  }
  x <- array(slope * 1000 / 52, dim = c(nrow(slope), ncol(slope), 1),
             dimnames = list(NULL, NULL, "standardized_slope"))
  as.data.frame(posterior::summarise_draws(posterior::as_draws_array(x),
    "rhat", "ess_bulk", "ess_tail", "mcse_mean", "sd"))
}

diagnostic_summary <- function(fit, prepared, model, variant, label) {
  d <- diagnostics(fit)
  ds <- derived_diagnostics(fit, prepared, model)
  p <- rbind(d$parameters[, names(ds)], ds)
  p <- p[!p$variable %in% c("L_Omega[1,1]", "L_Omega[1,2]"), ]
  ok <- d$passed && all(is.finite(ds$rhat)) && all(ds$rhat < 1.01) &&
    all(ds$ess_bulk >= 400) && all(ds$ess_tail >= 400)
  data.frame(run = label, model = model, variant = variant,
    max_rhat = max(p$rhat), min_bulk_ess = min(p$ess_bulk), min_tail_ess = min(p$ess_tail),
    max_slope_mcse = max(ds$mcse_mean), divergences = sum(d$sampler$divergences),
    treedepth_hits = sum(d$sampler$treedepth_hits), min_ebfmi = min(d$sampler$ebfmi), passed = ok)
}

population_curves <- function(draws, prepared, model) {
  weeks <- seq(min(prepared$obs$week), max(prepared$obs$week), length.out = 70)
  covs <- colMeans(prepared$patients[, c("age10", "male", "current", "never")])
  profiles <- list(cohort_standardized = covs)
  if (model == "C") profiles <- c(profiles, list(
    female_ex_smoker_age65 = c(0, 0, 0, 0), male_ex_smoker_age65 = c(0, 1, 0, 0),
    female_never_smoker_age65 = c(0, 0, 0, 1), female_current_smoker_age65 = c(0, 0, 1, 0)))
  do.call(rbind, lapply(names(profiles), function(name) {
    v <- unname(profiles[[name]])
    grid <- data.frame(time = weeks / 52, age10 = v[1], male = v[2], current = v[3], never = v[4])
    mu <- draws$beta %*% t(design_matrix(grid, model)) * 1000
    data.frame(model = model, time_week = weeks, profile = name,
      mean = colMeans(mu), lower95 = apply(mu, 2, quantile, .025), upper95 = apply(mu, 2, quantile, .975))
  }))
}

prior_predictive <- function(prepared, model, seed = 4224L, n_draws = 1000L) {
  set.seed(seed)
  X <- design_matrix(prepared$obs, model)
  scales <- if (model == "C") c(1.5, .5, .5, 1, 1, 1, .25, .5, .5, .5) else c(1.5, .5)
  beta <- matrix(rnorm(n_draws * ncol(X)), n_draws) * rep(scales, each = n_draws)
  beta[, 1] <- beta[, 1] + 3
  sigma <- abs(rnorm(n_draws, 0, .5))
  pred <- beta %*% t(X)
  if (model != "P") {
    tau1 <- abs(rnorm(n_draws, 0, 1)); tau2 <- abs(rnorm(n_draws, 0, .5))
    rho <- 2 * rbeta(n_draws, 2, 2) - 1
    for (j in seq_len(nrow(prepared$patients))) {
      z1 <- rnorm(n_draws); z2 <- rnorm(n_draws)
      i <- which(prepared$obs$id == j)
      pred[, i] <- pred[, i, drop = FALSE] + tau1 * z1 +
        outer(tau2 * (rho * z1 + sqrt(1 - rho^2) * z2), prepared$obs$time[i])
    }
  }
  pred <- pred + matrix(rnorm(length(pred)), nrow(pred)) * sigma
  stats <- cbind(negative_fraction = rowMeans(pred < 0),
    median_fvc_ml = apply(pred, 1, median) * 1000,
    q01_fvc_ml = apply(pred, 1, quantile, .01) * 1000,
    q99_fvc_ml = apply(pred, 1, quantile, .99) * 1000)
  obs <- c(mean(prepared$obs$y < 0), median(prepared$obs$y) * 1000,
    quantile(prepared$obs$y, .01) * 1000, quantile(prepared$obs$y, .99) * 1000)
  do.call(rbind, lapply(seq_len(ncol(stats)), function(i)
    data.frame(model = model, statistic = colnames(stats)[i], observed = unname(obs[i]),
      lower = unname(quantile(stats[, i], .025)), median = median(stats[, i]),
      upper = unname(quantile(stats[, i], .975)), ppp = mean(stats[, i] >= obs[i]))))
}

aggregate_cv <- function(cv) {
  metrics <- do.call(rbind, lapply(split(cv, cv$model), function(d)
    data.frame(model = d$model[1], mean_log_score = mean(d$log_score),
      rmse = sqrt(mean(d$mse)), mae = mean(d$mae), coverage95 = mean(d$coverage95), width95 = mean(d$width95))))
  comparisons <- do.call(rbind, lapply(list(c("H", "P"), c("C", "P"), c("C", "H")), function(pair) {
    a <- cv[cv$model == pair[1], ]; b <- cv[cv$model == pair[2], ]
    difference <- a$log_score - b$log_score[match(a$patient, b$patient)]
    data.frame(model_a = pair[1], model_b = pair[2], difference = mean(difference),
      se = sd(difference) / sqrt(length(difference)))
  }))
  list(metrics = metrics, comparisons = comparisons)
}
