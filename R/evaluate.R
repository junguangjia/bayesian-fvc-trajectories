# Predictive evaluation for previously unseen patients. All internal model
# quantities use liters and years; evaluation summaries use milliliters.

eval_log_mean_exp <- function(x) {
  if (!length(x) || anyNA(x)) stop("Log densities must be nonempty and nonmissing.")
  largest <- max(x)
  if (is.infinite(largest)) return(largest)
  largest + log(mean(exp(x - largest)))
}

eval_draw_indices <- function(n, maximum = 1000L) {
  if (n <= maximum) seq_len(n) else unique(round(seq(1, n, length.out = maximum)))
}

eval_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(seed)
  force(code)
}

eval_validate_draws <- function(draws, model, k) {
  if (!model %in% c("P", "H", "C")) stop("model must be P, H, or C.")
  if (!is.matrix(draws$beta) || ncol(draws$beta) != k ||
      !nrow(draws$beta) || any(!is.finite(draws$beta))) {
    stop("draws$beta must be a finite draws-by-coefficients matrix.")
  }
  s <- nrow(draws$beta)
  if (length(draws$sigma) != s || any(!is.finite(draws$sigma)) ||
      any(draws$sigma <= 0)) stop("Each posterior draw needs a positive sigma.")
  if (model != "P") eval_random_factor(draws, seq_len(s))
  invisible(s)
}

# The three columns encode the nonzero elements of diag(tau) %*% L_Omega.
# A rho vector may be supplied when the correlation Cholesky factor is absent.
eval_random_factor <- function(draws, indices) {
  s <- nrow(draws$beta)
  if (!is.matrix(draws$tau) || !identical(dim(draws$tau), c(s, 2L)) ||
      any(!is.finite(draws$tau)) || any(draws$tau < 0)) {
    stop("Hierarchical draws need a nonnegative draws-by-2 tau matrix.")
  }
  if (!is.null(draws$L_Omega)) {
    if (!identical(dim(draws$L_Omega), c(s, 2L, 2L)) ||
        any(!is.finite(draws$L_Omega))) stop("Invalid L_Omega dimensions or values.")
    l <- draws$L_Omega
    if (any(abs(l[, 1, 2]) > 1e-8) || any(abs(l[, 1, 1] - 1) > 1e-8) ||
        any(l[, 2, 2] < 0) ||
        any(abs(l[, 2, 1]^2 + l[, 2, 2]^2 - 1) > 1e-7)) {
      stop("L_Omega must be a lower Cholesky factor of a correlation matrix.")
    }
    correlation <- l[indices, 2, 1]
    diagonal <- l[indices, 2, 2]
  } else {
    if (length(draws$rho) != s || any(!is.finite(draws$rho)) ||
        any(abs(draws$rho) > 1)) stop("Hierarchical draws need L_Omega or rho.")
    correlation <- draws$rho[indices]
    diagonal <- sqrt(pmax(0, 1 - correlation^2))
  }
  cbind(b11 = draws$tau[indices, 1],
        b21 = draws$tau[indices, 2] * correlation,
        b22 = draws$tau[indices, 2] * diagonal)
}

# Integrate the two shared patient random effects analytically. The only
# inverse is 2-by-2, regardless of the number of visits. The final Jacobian
# converts a joint density in liters to the corresponding density in mL.
eval_joint_log_density <- function(y, time, mu, sigma, factor = NULL,
                                   units_scale = 1000) {
  if (!is.matrix(mu) || ncol(mu) != length(y) || length(time) != length(y) ||
      nrow(mu) != length(sigma) || any(!is.finite(mu)) ||
      any(!is.finite(y)) || any(!is.finite(time)) ||
      any(!is.finite(sigma)) || any(sigma <= 0) ||
      length(units_scale) != 1L || !is.finite(units_scale) || units_scale <= 0) {
    stop("Invalid Gaussian predictive-density inputs.")
  }
  n <- length(y)
  residual <- sweep(-mu, 2, y, "+")
  variance <- sigma^2
  quadratic <- rowSums(residual^2) / variance
  log_determinant <- n * log(variance)
  if (!is.null(factor)) {
    if (!is.matrix(factor) || !identical(dim(factor), c(nrow(mu), 3L)) ||
        any(!is.finite(factor))) stop("Invalid random-effect factor.")
    b11 <- factor[, 1]
    b21 <- factor[, 2]
    b22 <- factor[, 3]
    t1 <- sum(time)
    t2 <- sum(time^2)
    a11 <- 1 + (n * b11^2 + 2 * t1 * b11 * b21 + t2 * b21^2) / variance
    a12 <- (t1 * b11 * b22 + t2 * b21 * b22) / variance
    a22 <- 1 + t2 * b22^2 / variance
    # Expanded determinant avoids cancellation between the leading products.
    determinant <- 1 +
      (n * b11^2 + 2 * t1 * b11 * b21 + t2 * (b21^2 + b22^2)) / variance +
      b11^2 * b22^2 * pmax(0, n * t2 - t1^2) / variance^2
    r1 <- rowSums(residual)
    r2 <- as.vector(residual %*% time)
    w1 <- b11 * r1 + b21 * r2
    w2 <- b22 * r2
    correction <- (a22 * w1^2 - 2 * a12 * w1 * w2 + a11 * w2^2) /
      (determinant * variance^2)
    quadratic <- pmax(0, quadratic - correction)
    log_determinant <- log_determinant + log(determinant)
  }
  -0.5 * (n * log(2 * pi) + log_determinant + quadratic) - n * log(units_scale)
}

evaluate_new_patients <- function(draws, prepared, model, seed = 4224) {
  obs <- prepared$obs
  required <- c("patient", "y", "time")
  if (!all(required %in% names(obs)) || !nrow(obs) || anyNA(obs$patient)) {
    stop("prepared$obs must contain nonempty patient, y, and time columns.")
  }
  x <- design_matrix(obs, model)
  s <- eval_validate_draws(draws, model, ncol(x))
  selected <- eval_draw_indices(s)
  m <- length(selected)
  random_factor <- if (model == "P") NULL else eval_random_factor(draws, seq_len(s))
  groups <- split(seq_len(nrow(obs)), as.character(obs$patient))
  rows <- eval_with_seed(seed, lapply(groups, function(index) {
    local_x <- x[index, , drop = FALSE]
    mu <- draws$beta %*% t(local_x)
    log_density <- eval_joint_log_density(
      obs$y[index], obs$time[index], mu, draws$sigma, random_factor
    )
    prediction <- mu[selected, , drop = FALSE]
    if (model != "P") {
      factor <- random_factor[selected, , drop = FALSE]
      z1 <- rnorm(m)
      z2 <- rnorm(m)
      intercept <- factor[, 1] * z1
      slope <- factor[, 2] * z1 + factor[, 3] * z2
      prediction <- prediction + intercept + slope %o% obs$time[index]
    }
    prediction <- prediction + matrix(rnorm(m * length(index)), m) * draws$sigma[selected]
    intervals <- apply(prediction * 1000, 2, quantile, probs = c(0.025, 0.975))
    actual <- obs$y[index] * 1000
    error <- (colMeans(mu) - obs$y[index]) * 1000
    data.frame(
      patient = as.character(obs$patient[index[1]]),
      log_score = eval_log_mean_exp(log_density),
      mse = mean(error^2), mae = mean(abs(error)),
      coverage95 = mean(actual >= intervals[1, ] & actual <= intervals[2, ]),
      width95 = mean(intervals[2, ] - intervals[1, ]),
      stringsAsFactors = FALSE
    )
  }))
  per_patient <- do.call(rbind, rows)
  rownames(per_patient) <- NULL
  list(per_patient = per_patient, metadata = list(
    total_draws = s, predictive_draws = m, seed = seed,
    density_units = "mL", squared_error_units = "mL^2",
    interval_type = "pointwise 95% posterior predictive intervals",
    patient_effects = "integrated for scores; one shared effect per patient per simulated draw"
  ))
}

eval_safe_sd <- function(x) if (length(x) > 1L) sd(x) else NA_real_
eval_safe_cor <- function(x, y) {
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) NA_real_ else cor(x, y)
}

eval_check_context <- function(obs) {
  groups <- split(seq_len(nrow(obs)), as.character(obs$patient))
  means <- matrix(0, nrow(obs), length(groups))
  slopes <- matrix(0, nrow(obs), length(groups))
  valid_slope <- logical(length(groups))
  previous <- following <- integer()
  for (j in seq_along(groups)) {
    ix <- groups[[j]]
    means[ix, j] <- 1 / length(ix)
    centered <- obs$time[ix] - mean(obs$time[ix])
    denominator <- sum(centered^2)
    if (denominator > 0) {
      slopes[ix, j] <- centered / denominator
      valid_slope[j] <- TRUE
    }
    ordered <- ix[order(obs$time[ix])]
    if (length(ordered) > 1L) {
      previous <- c(previous, head(ordered, -1L))
      following <- c(following, tail(ordered, -1L))
    }
  }
  time_centered <- obs$time - mean(obs$time)
  time_denominator <- sum(time_centered^2)
  list(mean_weights = means, slope_weights = slopes[, valid_slope, drop = FALSE],
       previous = previous, following = following,
       time_weights = if (time_denominator > 0) time_centered / time_denominator else NULL)
}

eval_distribution_statistics <- function(y, context) {
  patient_means <- as.vector(y %*% context$mean_weights)
  patient_slopes <- as.vector(y %*% context$slope_weights)
  c(median_mL = median(y), q05_mL = unname(quantile(y, 0.05)),
    q95_mL = unname(quantile(y, 0.95)), sd_mL = eval_safe_sd(y),
    patient_mean_sd_mL = eval_safe_sd(patient_means),
    patient_slope_mean_mL_per_year = if (length(patient_slopes)) mean(patient_slopes) else NA_real_,
    patient_slope_sd_mL_per_year = eval_safe_sd(patient_slopes))
}

eval_residual_statistics <- function(residual, context) {
  c(residual_lag1_correlation = eval_safe_cor(
      residual[context$previous], residual[context$following]),
    residual_time_slope_mL_per_year = if (is.null(context$time_weights)) NA_real_
      else sum(residual * context$time_weights),
    absolute_residual_time_slope_mL_per_year = if (is.null(context$time_weights)) NA_real_
      else sum(abs(residual) * context$time_weights))
}

posterior_checks <- function(draws, prepared, model, seed = 4224) {
  obs <- prepared$obs
  x <- design_matrix(obs, model)
  s <- eval_validate_draws(draws, model, ncol(x))
  selected <- eval_draw_indices(s)
  m <- length(selected)
  mu <- draws$beta[selected, , drop = FALSE] %*% t(x)
  if (model != "P") {
    dimensions <- dim(draws$u)
    if (length(dimensions) != 3L || dimensions[1] != s || dimensions[2] != 2L ||
        any(!is.finite(draws$u)) || anyNA(obs$id) || any(obs$id < 1) ||
        any(obs$id != as.integer(obs$id)) || any(obs$id > dimensions[3])) {
      stop("Conditional posterior checks need patient effects u[draw, effect, patient].")
    }
    intercept <- matrix(draws$u[selected, 1, obs$id, drop = FALSE], m, nrow(obs))
    slope <- matrix(draws$u[selected, 2, obs$id, drop = FALSE], m, nrow(obs))
    mu <- mu + intercept + sweep(slope, 2, obs$time, "*")
  }
  replicated <- eval_with_seed(seed, mu + matrix(rnorm(m * nrow(obs)), m) * draws$sigma[selected])
  context <- eval_check_context(obs)
  actual <- obs$y * 1000
  observed_distribution <- eval_distribution_statistics(actual, context)
  observed_rows <- replicated_rows <- vector("list", m)
  for (i in seq_len(m)) {
    fitted <- mu[i, ] * 1000
    predicted <- replicated[i, ] * 1000
    observed_rows[[i]] <- c(observed_distribution,
      eval_residual_statistics(actual - fitted, context))
    replicated_rows[[i]] <- c(eval_distribution_statistics(predicted, context),
      eval_residual_statistics(predicted - fitted, context))
  }
  observed <- do.call(rbind, observed_rows)
  replicated_statistics <- do.call(rbind, replicated_rows)
  result <- lapply(seq_len(ncol(observed)), function(k) {
    valid <- is.finite(observed[, k]) & is.finite(replicated_statistics[, k])
    if (any(valid)) {
      limits <- quantile(replicated_statistics[valid, k], c(0.025, 0.5, 0.975), names = FALSE)
      observed_value <- median(observed[valid, k])
      ppp <- mean(replicated_statistics[valid, k] >= observed[valid, k])
    } else {
      limits <- rep(NA_real_, 3)
      observed_value <- ppp <- NA_real_
    }
    data.frame(statistic = colnames(observed)[k], observed = observed_value,
               lower = limits[1], median = limits[2], upper = limits[3], ppp = ppp)
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  attr(result, "metadata") <- list(
    total_draws = s, predictive_draws = m, seed = seed,
    observation_units = "mL", time_units = "years",
    conditioning = "observed patients and their posterior random effects",
    observed_residual_statistics = "median of draw-specific discrepancies",
    ppp = "Pr(replicated discrepancy >= observed discrepancy), paired within each draw"
  )
  result
}
