# Run from the project root with Rscript --vanilla tests/test-model.R.
# These tests exercise pure input and diagnostic logic without compiling Stan.
source("R/model.R")

design_matrix <- function(obs, model) {
  x <- cbind(alpha = 1, beta_time = obs$time)
  if (model == "C") {
    x <- cbind(x, beta_age = obs$age10, beta_male = obs$male,
               beta_current = obs$current, beta_never = obs$never,
               gamma_age = obs$time * obs$age10,
               gamma_male = obs$time * obs$male,
               gamma_current = obs$time * obs$current,
               gamma_never = obs$time * obs$never)
  }
  x
}

expect_error <- function(expression) {
  failed <- tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
  stopifnot(failed)
}

obs <- data.frame(id = rep(1:3, each = 3),
                  patient = rep(c("synthetic-a", "synthetic-b", "synthetic-c"), each = 3),
                  y = c(3, 2.9, 2.8, 2.5, 2.4, 2.3, 4, 3.9, 3.8),
                  time = rep(c(0, 0.5, 1), 3),
                  age10 = rep(c(-1, 0, 1), each = 3),
                  male = rep(c(0, 1, 1), each = 3),
                  current = rep(c(0, 0, 1), each = 3),
                  never = rep(c(0, 1, 0), each = 3))
prepared <- list(obs = obs, patients = obs[!duplicated(obs$id), ],
                 audit = list(synthetic = TRUE))
primary <- .fvc_model_inputs(prepared, "C", 1)
wide <- .fvc_model_inputs(prepared, "C", 2)
stopifnot(primary$data$K == 10L, primary$data$J == 3L,
          identical(primary$data$patient, obs$id),
          identical(primary$data$beta_prior_mean, c(3, rep(0, 9))),
          identical(wide$data$beta_prior_sd, 2 * primary$data$beta_prior_sd),
          identical(wide$data$tau_prior_sd, 2 * primary$data$tau_prior_sd),
          identical(wide$data$sigma_prior_sd, 2 * primary$data$sigma_prior_sd))
stopifnot(is.null(.fvc_model_inputs(prepared, "P", 1)$data$J),
          .fvc_model_inputs(prepared, "H", 1)$data$K == 2L)
broken <- prepared
broken$obs$id[1] <- 4L
expect_error(.fvc_model_inputs(broken, "H", 1))
broken <- prepared
broken$obs$y[1] <- NA_real_
expect_error(.fvc_model_inputs(broken, "P", 1))
expect_error(.fvc_settings(list(iter = 100, warmup = 100)))
expect_error(.fvc_settings(list(adapt_delta = 1)))
expect_error(.fvc_settings(list(unrecognized = TRUE)))
stopifnot(identical(.fvc_settings(list())$save_warmup, FALSE),
          .fvc_settings(list(chains = 2))$cores == 2L)

# Exercise the production cache key, not a separately reconstructed test hash.
cache_settings <- .fvc_settings(list())
cache_software <- c(R = "test-R", rstan = "test-rstan", RcppParallel = "test-rpp",
                    QuickJSR = "test-qjs")
key <- function(p = prepared, input = primary, setting = cache_settings,
                 seed = 4224L, source = "stan-source-a", implementation = "R-source-a",
                 software = cache_software) {
  .fvc_cache_identity(p, input, setting, seed, source, implementation, software)
}
baseline <- key()
changed_data <- prepared
changed_data$obs$y[1] <- changed_data$obs$y[1] + 0.001
changed_ids <- prepared
changed_ids$obs$patient[1] <- "synthetic-different-identifier"
changed_setting <- cache_settings
changed_setting$adapt_delta <- 0.99
changed_software <- cache_software
changed_software["QuickJSR"] <- "test-qjs-updated"
different <- list(key(p = changed_data), key(p = changed_ids), key(input = wide),
                  key(setting = changed_setting), key(seed = 4225L),
                  key(source = "stan-source-b"), key(implementation = "R-source-b"),
                  key(software = changed_software))
stopifnot(identical(key(), baseline),
          all(vapply(different, function(x) !identical(x$content_sha256,
                                                       baseline$content_sha256), logical(1))),
          !identical(key(source = "stan-source-b")$compiled_sha256, baseline$compiled_sha256),
          !identical(key(implementation = "R-source-b")$compiled_sha256, baseline$compiled_sha256),
          !identical(key(software = changed_software)$compiled_sha256, baseline$compiled_sha256),
          all(c("RcppParallel", "QuickJSR") %in% names(.fvc_software())))
# The implementation identity is captured at source time rather than rereading
# its pathname while a long research run is already in progress.
old_path <- .fvc_model_file
.fvc_model_file <- "deliberately-nonexistent-model-source.R"
captured <- .fvc_cache_identity(prepared, primary, cache_settings, 4224L,
                               "stan-source-a", software = cache_software)
.fvc_model_file <- old_path
stopifnot(identical(captured$identity$implementation_sha256,
                    .fvc_model_implementation_sha256))

normalized <- .fvc_normalize_parameter_draws(list(
  sigma = array(1:3, 3), rho = array(c(-.1, 0, .1), 3), lp__ = array(-1:-3, 3),
  beta = matrix(1:6, 3), tau = matrix(1:6, 3)))
stopifnot(is.null(dim(normalized$sigma)), is.null(dim(normalized$rho)),
          is.null(dim(normalized$lp__)), identical(dim(normalized$beta), c(3L, 2L)),
          identical(dim(normalized$tau), c(3L, 2L)))
# Session-cache lookup retains the same object even when no disk file exists.
sentinel <- new.env(parent = emptyenv())
assign("test-retained-reference", list(model = sentinel, module = sentinel),
       envir = .fvc_compiled_models)
gc(verbose = FALSE)
stopifnot(identical(.fvc_get_compiled_model("unused.stan", "test-retained-reference",
                                           "unused-cache"), sentinel))
rm(list = "test-retained-reference", envir = .fvc_compiled_models)

set.seed(1101)
draws <- array(rnorm(2000 * 4 * 3), dim = c(2000, 4, 3),
               dimnames = list(NULL, NULL, c("beta[1]", "sigma", "L_Omega[1,1]")))
draws[, , "sigma"] <- exp(draws[, , "sigma"])
draws[, , "L_Omega[1,1]"] <- 1
good <- .fvc_parameter_diagnostics(draws)
stopifnot(good$passed, identical(names(good$parameters),
                                c("variable", "rhat", "ess_bulk", "ess_tail", "mcse_mean", "sd")))
separated <- draws
separated[, 1, "beta[1]"] <- separated[, 1, "beta[1]"] + 5
stopifnot(!.fvc_parameter_diagnostics(separated)$passed)
stuck <- draws
stuck[, , "beta[1]"] <- 1
stopifnot(!.fvc_parameter_diagnostics(stuck)$passed)
nonfinite <- draws
nonfinite[1, 1, "beta[1]"] <- NA_real_
stopifnot(!.fvc_parameter_diagnostics(nonfinite)$passed)

samplers <- lapply(1:4, function(i) {
  cbind(divergent__ = rep(0, 2000), treedepth__ = rep(4, 2000),
        energy__ = rnorm(2000))
})
stopifnot(.fvc_sampler_diagnostics(samplers, rep(12L, 4))$passed)
divergent <- samplers
divergent[[1]][1, "divergent__"] <- 1
stopifnot(!.fvc_sampler_diagnostics(divergent, rep(12L, 4))$passed)
saturated <- samplers
saturated[[1]][1, "treedepth__"] <- 12
stopifnot(!.fvc_sampler_diagnostics(saturated, rep(12L, 4))$passed)
low_energy <- samplers
low_energy[[1]][, "energy__"] <- seq_len(2000)
stopifnot(!.fvc_sampler_diagnostics(low_energy, rep(12L, 4))$passed)
cat("Model inputs, cache invalidation, scalar normalization, and diagnostic failure gates passed.\n")
