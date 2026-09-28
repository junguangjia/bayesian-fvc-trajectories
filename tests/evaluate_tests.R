run_evaluate_tests <- function() {
  close <- function(actual, expected, tolerance = 1e-9) {
    stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected),
                               tolerance = tolerance)))
  }
  n <- 7L
  s <- 9L
  time <- seq(-0.2, 1.3, length.out = n)
  y <- c(2.8, 2.9, 2.6, 2.7, 2.5, 2.6, 2.3)
  mu <- outer(seq(2.7, 3.1, length.out = s), rep(1, n)) -
    outer(seq(0.1, 0.4, length.out = s), time)
  sigma <- seq(0.08, 0.25, length.out = s)
  tau <- cbind(seq(0, 0.7, length.out = s), seq(0.4, 0, length.out = s))
  rho <- seq(-0.95, 0.95, length.out = s)
  draws <- list(beta = cbind(seq(2.7, 3.1, length.out = s), -0.2),
                sigma = sigma, tau = tau, rho = rho)
  factor <- eval_random_factor(draws, seq_len(s))
  actual <- eval_joint_log_density(y, time, mu, sigma, factor)
  z <- cbind(1, time)
  expected <- vapply(seq_len(s), function(i) {
    l <- matrix(c(factor[i, 1], factor[i, 2], 0, factor[i, 3]), 2)
    covariance <- diag(sigma[i]^2, n) + z %*% tcrossprod(l) %*% t(z)
    residual <- y - mu[i, ]
    -0.5 * (n * log(2 * pi) + as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
            as.numeric(crossprod(residual, solve(covariance, residual)))) - n * log(1000)
  }, numeric(1))
  close(actual, expected)

  # Exactly zero random-effect scales recover the independent Gaussian model.
  pooled <- eval_joint_log_density(y, time, mu, sigma)
  zero_random <- eval_joint_log_density(y, time, mu, sigma, matrix(0, s, 3))
  close(pooled, zero_random)
  expected_pooled <- vapply(seq_len(s), function(i) {
    sum(dnorm(y * 1000, mu[i, ] * 1000, sigma[i] * 1000, log = TRUE))
  }, numeric(1))
  close(pooled, expected_pooled)
  close(eval_joint_log_density(y, time, mu, sigma, factor, units_scale = 1) - actual,
        rep(n * log(1000), s))
  close(eval_log_mean_exp(c(-1000000, -1000001)),
        -1000000 + log((1 + exp(-1)) / 2))
  stopifnot(identical(eval_log_mean_exp(c(-Inf, -Inf)), -Inf))

  # A supplied Cholesky factor and the corresponding rho give the same factors.
  draws$L_Omega <- array(0, c(s, 2L, 2L))
  draws$L_Omega[, 1, 1] <- 1
  draws$L_Omega[, 2, 1] <- rho
  draws$L_Omega[, 2, 2] <- sqrt(1 - rho^2)
  close(eval_random_factor(draws, seq_len(s)), factor)

  # Bind only the design function locally, leaving the application's function intact.
  test_environment <- new.env(parent = environment(evaluate_new_patients))
  test_environment$design_matrix <- function(obs, model) cbind(intercept = 1, time = obs$time)
  evaluate <- evaluate_new_patients
  checks <- posterior_checks
  environment(evaluate) <- test_environment
  environment(checks) <- test_environment
  test_obs <- data.frame(
    id = rep(1:3, each = 4), patient = rep(c("synthetic_a", "synthetic_b", "synthetic_c"), each = 4),
    y = c(3.1, 2.9, 3.0, 2.7, 2.7, 2.6, 2.8, 2.4, 3.4, 3.2, 3.1, 3.0),
    time = rep(c(0, 0.3, 0.6, 1), 3)
  )
  prepared <- list(obs = test_obs)
  many <- 1201L
  test_draws <- list(beta = cbind(rep(3, many), rep(-0.25, many)),
                     sigma = rep(0.15, many), tau = cbind(rep(0.35, many), rep(0.1, many)),
                     rho = rep(0.3, many), u = array(0, c(many, 2L, 3L)))
  set.seed(773)
  rng_before <- .Random.seed
  result <- evaluate(test_draws, prepared, "H", seed = 193)
  stopifnot(identical(.Random.seed, rng_before), nrow(result$per_patient) == 3L,
            result$metadata$total_draws == many, result$metadata$predictive_draws == 1000L,
            all(is.finite(as.matrix(result$per_patient[, -1]))),
            identical(result, evaluate(test_draws, prepared, "H", seed = 193)))
  altered_effects <- test_draws
  altered_effects$u[] <- 100
  # New-patient evaluation must never reuse fitted effects for observed patients.
  stopifnot(identical(result, evaluate(altered_effects, prepared, "H", seed = 193)))
  expected_mse <- mean(((3 - 0.25 * test_obs$time[1:4]) - test_obs$y[1:4])^2) * 1e6
  close(result$per_patient$mse[1], expected_mse)

  check_result <- checks(test_draws, prepared, "H", seed = 193)
  stopifnot(identical(names(check_result), c("statistic", "observed", "lower", "median", "upper", "ppp")),
            nrow(check_result) == 10L, all(is.finite(as.matrix(check_result[, -1]))),
            all(check_result$ppp >= 0 & check_result$ppp <= 1),
            !any(test_obs$patient %in% unlist(check_result)),
            identical(.Random.seed, rng_before))
  context <- eval_check_context(test_obs)
  stopifnot(length(context$previous) == 9L, length(context$following) == 9L,
            all(test_obs$patient[context$previous] == test_obs$patient[context$following]))

  # Single-visit patients have well-defined scores and intervals; slope checks are NA.
  single <- list(obs = test_obs[c(1, 5, 9), ])
  single_result <- evaluate(test_draws, single, "P", seed = 193)
  stopifnot(nrow(single_result$per_patient) == 3L,
            all(is.finite(as.matrix(single_result$per_patient[, -1]))))
  single_checks <- checks(test_draws, single, "P", seed = 193)
  slope_rows <- grepl("patient_slope_", single_checks$statistic)
  stopifnot(all(is.na(single_checks$observed[slope_rows])))
  invisible(TRUE)
}
