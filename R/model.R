# Model fitting and diagnostics. Outcomes are liters and time is years.
# Source this file after defining design_matrix(obs, model).

.fvc_model_file <- local({
  candidates <- lapply(sys.frames(), function(frame) frame$ofile)
  candidates <- Filter(function(x) is.character(x) && length(x) == 1L,
                       candidates)
  path <- if (length(candidates)) tail(candidates, 1L)[[1L]] else "R/model.R"
  normalizePath(path, mustWork = TRUE)
})

# Capture the implementation that was actually sourced. Later edits on disk must
# not relabel functions already running in this R session with a different hash.
.fvc_model_implementation_sha256 <- digest::digest(file = .fvc_model_file,
                                                  algo = "sha256")
.fvc_compiled_models <- new.env(parent = emptyenv())

.fvc_require_packages <- function() {
  required <- c("rstan", "posterior", "digest")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (!all(available)) {
    stop("Missing required packages: ", paste(required[!available], collapse = ", "),
         call. = FALSE)
  }
}

.fvc_scalar <- function(x, name, lower, upper = Inf, integer = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= lower && x <= upper && (!integer || x == as.integer(x))
  if (!valid) stop("Invalid setting: ", name, call. = FALSE)
  if (integer) as.integer(x) else as.numeric(x)
}

.fvc_settings <- function(settings) {
  defaults <- list(chains = 4L, iter = 6000L, warmup = 3000L,
                   adapt_delta = 0.95, max_treedepth = 12L, cores = 4L,
                   refresh = 0L)
  if (!is.list(settings) || (length(settings) && is.null(names(settings)))) {
    stop("settings must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(settings), c(names(defaults), "stan_dir"))
  if (length(unknown)) stop("Unknown settings: ", paste(unknown, collapse = ", "),
                            call. = FALSE)
  for (name in intersect(names(settings), names(defaults))) {
    defaults[[name]] <- settings[[name]]
  }
  for (name in c("chains", "iter", "warmup", "max_treedepth", "cores")) {
    defaults[[name]] <- .fvc_scalar(defaults[[name]], name, 1, integer = TRUE)
  }
  defaults$refresh <- .fvc_scalar(defaults$refresh, "refresh", 0, integer = TRUE)
  defaults$adapt_delta <- .fvc_scalar(defaults$adapt_delta, "adapt_delta", 0, 1)
  if (defaults$adapt_delta <= 0 || defaults$adapt_delta >= 1) {
    stop("adapt_delta must be strictly between zero and one.", call. = FALSE)
  }
  if (defaults$iter <= defaults$warmup) {
    stop("iter must exceed warmup.", call. = FALSE)
  }
  defaults$cores <- min(defaults$cores, defaults$chains)
  defaults$thin <- 1L
  defaults$save_warmup <- FALSE
  defaults
}

.fvc_software <- function() {
  packages <- c("rstan", "StanHeaders", "Rcpp", "RcppEigen", "RcppParallel",
                "QuickJSR", "BH", "Matrix",
                "posterior", "digest")
  versions <- vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      "unavailable"
    }
  }, character(1))
  c(R = as.character(getRversion()), platform = R.version$platform, versions)
}

.fvc_model_inputs <- function(prepared, model, prior_scale) {
  model <- match.arg(model, c("P", "H", "C"))
  prior_scale <- .fvc_scalar(prior_scale, "prior_scale", .Machine$double.eps)
  if (!is.list(prepared) || !is.data.frame(prepared$obs) ||
      !is.data.frame(prepared$patients)) {
    stop("prepared must contain obs and patients data frames.", call. = FALSE)
  }
  obs <- prepared$obs
  required <- c("id", "patient", "y", "time", "age10", "male", "current", "never")
  if (!all(required %in% names(obs)) || !nrow(obs)) {
    stop("Observation data do not have the required nonempty schema.", call. = FALSE)
  }
  numeric_columns <- setdiff(required, "patient")
  if (!all(vapply(obs[numeric_columns], is.numeric, logical(1))) ||
      !all(is.finite(as.matrix(obs[numeric_columns])))) {
    stop("All numeric observation fields must be finite.", call. = FALSE)
  }
  J <- nrow(prepared$patients)
  if (!J || anyNA(obs$patient) || any(!nzchar(as.character(obs$patient))) ||
      any(obs$id != as.integer(obs$id)) ||
      !identical(sort(unique(as.integer(obs$id))), seq_len(J))) {
    stop("Participant indices must cover 1:J with nonmissing private identifiers.",
         call. = FALSE)
  }
  mapping <- unique(obs[c("id", "patient")])
  if (nrow(mapping) != J || anyDuplicated(mapping$id) ||
      anyDuplicated(mapping$patient)) {
    stop("Each participant index must map to exactly one participant.", call. = FALSE)
  }
  if (any(!as.matrix(obs[c("male", "current", "never")]) %in% c(0, 1)) ||
      any(obs$current + obs$never > 1)) {
    stop("Sex and smoking indicators must use the documented binary encoding.",
         call. = FALSE)
  }
  if (!exists("design_matrix", mode = "function", inherits = TRUE)) {
    stop("Source the preprocessing module defining design_matrix first.", call. = FALSE)
  }
  X <- design_matrix(obs, model)
  all_names <- c("alpha", "beta_time", "beta_age", "beta_male", "beta_current",
                 "beta_never", "gamma_age", "gamma_male", "gamma_current",
                 "gamma_never")
  expected_names <- if (model == "C") all_names else all_names[1:2]
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) != nrow(obs) ||
      !identical(colnames(X), expected_names) || !all(is.finite(X))) {
    stop("design_matrix returned an invalid matrix or unexpected column order.",
         call. = FALSE)
  }
  if (any(X[, 1L] != 1) || !isTRUE(all.equal(unname(X[, 2L]), unname(obs$time)))) {
    stop("The first design columns must be the intercept and time in years.",
         call. = FALSE)
  }
  K <- ncol(X)
  prior_sd <- c(1.5, 0.5, 0.5, 1, 1, 1, 0.25, 0.5, 0.5, 0.5)[seq_len(K)] *
    prior_scale
  stan_data <- list(N = nrow(obs), K = K, X = unname(X), y = as.numeric(obs$y),
                    beta_prior_mean = c(3, rep(0, K - 1L)),
                    beta_prior_sd = prior_sd, sigma_prior_sd = 0.5 * prior_scale)
  if (model != "P") {
    stan_data$J <- J
    stan_data$patient <- as.integer(obs$id)
    stan_data$time <- as.numeric(obs$time)
    stan_data$tau_prior_sd <- c(1, 0.5) * prior_scale
  }
  list(model = model, data = stan_data, coefficient_names = expected_names,
       prior_scale = prior_scale, J = J)
}

.fvc_atomic_save <- function(object, path) {
  temporary <- tempfile(pattern = ".incomplete-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = "gzip")
  Sys.chmod(temporary, mode = "0600")
  if (!file.rename(temporary, path)) {
    stop("Could not atomically save the private cache file.", call. = FALSE)
  }
  invisible(path)
}

.fvc_cache_identity <- function(prepared, inputs, actual_settings, seed,
                                source_sha256,
                                implementation_sha256 = .fvc_model_implementation_sha256,
                                software = .fvc_software()) {
  data_hash <- digest::digest(list(obs = prepared$obs, patients = prepared$patients,
                                  audit = prepared$audit, stan_data = inputs$data),
                              algo = "sha256")
  settings_hash <- digest::digest(list(settings = actual_settings, seed = seed),
                                  algo = "sha256")
  software_hash <- digest::digest(software, algo = "sha256")
  identity <- list(version = 2L, model = inputs$model,
                   prior_scale = inputs$prior_scale, source_sha256 = source_sha256,
                   implementation_sha256 = implementation_sha256,
                   data_sha256 = data_hash, settings_sha256 = settings_hash,
                   software_sha256 = software_hash)
  list(identity = identity, content_sha256 = digest::digest(identity, algo = "sha256"),
       compiled_sha256 = digest::digest(list(source_sha256 = source_sha256,
                                             software_sha256 = software_hash,
                                             implementation_sha256 = implementation_sha256),
                                         algo = "sha256"))
}

.fvc_initialize_compiled <- function(model) {
  if (!inherits(model, "stanmodel")) stop("Invalid compiled-model object.")
  # mk_cppmodule triggers RStan's DSO loader and resolves its C++ class. A
  # deserialized object can have the right S4 class but invalid native pointers.
  module <- model@mk_cppmodule(model)
  if (is.null(module)) stop("Compiled Stan module could not be initialized.")
  list(model = model, module = module)
}

.fvc_get_compiled_model <- function(stan_file, compiled_hash, cache_dir) {
  if (exists(compiled_hash, envir = .fvc_compiled_models, inherits = FALSE)) {
    return(get(compiled_hash, envir = .fvc_compiled_models, inherits = FALSE)$model)
  }
  compiled_cache <- file.path(cache_dir, paste0("compiled-", compiled_hash, ".rds"))
  entry <- NULL
  if (file.exists(compiled_cache)) {
    entry <- tryCatch(.fvc_initialize_compiled(readRDS(compiled_cache)),
                      error = function(e) NULL)
    if (is.null(entry)) {
      message("Serialized compiled model could not initialize; recompiling current Stan source.")
    }
  }
  if (is.null(entry)) {
    # Supplying stanc_ret bypasses RStan's separate global/temp serialized-cache
    # lookup, which could otherwise return the same invalid object on recovery.
    parsed <- rstan::stanc(file = stan_file)
    if (!isTRUE(parsed$status)) stop("Stan source translation failed.", call. = FALSE)
    model <- rstan::stan_model(stanc_ret = parsed, auto_write = FALSE, save_dso = TRUE)
    entry <- .fvc_initialize_compiled(model)
    .fvc_atomic_save(model, compiled_cache)
  }
  # Retain both the model and resolved module even after fit objects are removed.
  assign(compiled_hash, entry, envir = .fvc_compiled_models)
  entry$model
}

fit_model <- function(prepared, model, prior_scale = 1, settings = list(),
                      cache_dir, seed = 4224L) {
  .fvc_require_packages()
  inputs <- .fvc_model_inputs(prepared, model, prior_scale)
  actual_settings <- .fvc_settings(settings)
  seed <- .fvc_scalar(seed, "seed", 1, .Machine$integer.max, integer = TRUE)
  stan_dir <- if (is.null(settings$stan_dir)) {
    file.path(dirname(dirname(.fvc_model_file)), "stan")
  } else settings$stan_dir
  source_name <- if (inputs$model == "P") "pooled.stan" else "hierarchical.stan"
  stan_file <- normalizePath(file.path(stan_dir, source_name), mustWork = TRUE)
  source_hash <- digest::digest(file = stan_file, algo = "sha256")
  software <- .fvc_software()
  cache_identity <- .fvc_cache_identity(prepared, inputs, actual_settings, seed,
                                       source_hash, software = software)
  identity <- cache_identity$identity
  content_hash <- cache_identity$content_sha256
  if (!is.character(cache_dir) || length(cache_dir) != 1L || !nzchar(cache_dir)) {
    stop("cache_dir must identify a private local cache directory.", call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  if (!dir.exists(cache_dir)) stop("Could not create cache directory.", call. = FALSE)
  fit_cache <- file.path(cache_dir, paste0("fit-", inputs$model, "-", content_hash, ".rds"))
  if (file.exists(fit_cache)) {
    cached <- readRDS(fit_cache)
    if (!is.list(cached) || !inherits(cached$fit, "stanfit") ||
        !identical(cached$metadata$content_sha256, content_hash)) {
      stop("Private fit cache failed its identity check.", call. = FALSE)
    }
    cached$metadata$cache_hit <- TRUE
    return(cached)
  }
  start <- proc.time()[["elapsed"]]
  stan_model <- .fvc_get_compiled_model(stan_file, cache_identity$compiled_sha256,
                                       cache_dir)
  fit <- rstan::sampling(
    stan_model, data = inputs$data, seed = seed,
    chains = actual_settings$chains, iter = actual_settings$iter,
    warmup = actual_settings$warmup, thin = actual_settings$thin,
    cores = actual_settings$cores, refresh = actual_settings$refresh,
    save_warmup = actual_settings$save_warmup,
    control = list(adapt_delta = actual_settings$adapt_delta,
                   max_treedepth = actual_settings$max_treedepth)
  )
  if (fit@mode != 0L) stop("Stan did not return a usable sampled fit.", call. = FALSE)
  metadata <- c(identity, list(
    content_sha256 = content_hash, source_file = source_name, software = software,
    settings = actual_settings, seed = seed, coefficient_names = inputs$coefficient_names,
    observations = inputs$data$N, participants = inputs$J,
    elapsed_seconds = unname(proc.time()[["elapsed"]] - start),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), cache_hit = FALSE
  ))
  result <- list(fit = fit, metadata = metadata)
  .fvc_atomic_save(result, fit_cache)
  result
}

extract_parameter_draws <- function(fit) {
  .fvc_normalize_parameter_draws(rstan::extract(fit, permuted = TRUE,
                                              inc_warmup = FALSE))
}

.fvc_normalize_parameter_draws <- function(draws) {
  # RStan may return scalar parameters as one-dimensional arrays. Explicit
  # vectors are needed for vector-by-matrix recycling in predictive operations.
  for (name in intersect(c("sigma", "rho", "lp__"), names(draws))) {
    draws[[name]] <- as.numeric(draws[[name]])
  }
  draws
}

.fvc_parameter_diagnostics <- function(draw_array) {
  if (length(dim(draw_array)) != 3L || is.null(dimnames(draw_array)[[3L]])) {
    stop("Expected an iteration-by-chain-by-parameter array.", call. = FALSE)
  }
  finite <- vapply(seq_len(dim(draw_array)[3L]), function(i) {
    all(is.finite(draw_array[, , i]))
  }, logical(1))
  constant <- vapply(seq_len(dim(draw_array)[3L]), function(i) {
    values <- as.vector(draw_array[, , i])
    all(is.finite(values)) && all(values == values[1L])
  }, logical(1))
  draws <- posterior::as_draws_array(draw_array)
  parameters <- as.data.frame(posterior::summarise_draws(
    draws, "rhat", "ess_bulk", "ess_tail", "mcse_mean", "sd"
  ))
  names(finite) <- names(constant) <- dimnames(draw_array)[[3L]]
  # Only deterministic Cholesky entries may legitimately remain constant.
  structural <- parameters$variable %in% c("L_Omega[1,1]", "L_Omega[1,2]") &
    constant[parameters$variable]
  checked <- !structural
  metrics <- as.matrix(parameters[checked, c("rhat", "ess_bulk", "ess_tail",
                                              "mcse_mean", "sd"), drop = FALSE])
  passed <- all(finite) && any(checked) && all(is.finite(metrics)) &&
    all(parameters$rhat[checked] < 1.01) &&
    all(parameters$ess_bulk[checked] >= 400) &&
    all(parameters$ess_tail[checked] >= 400) &&
    !any(constant[parameters$variable[checked]])
  list(parameters = parameters, passed = isTRUE(passed))
}

.fvc_sampler_diagnostics <- function(samplers, treedepth_limits) {
  if (length(samplers) != length(treedepth_limits) || !length(samplers)) {
    stop("Sampler output and per-chain treedepth limits must align.", call. = FALSE)
  }
  rows <- lapply(seq_along(samplers), function(chain) {
    samples <- samplers[[chain]]
    required <- c("divergent__", "treedepth__", "energy__")
    if (!is.matrix(samples) || !all(required %in% colnames(samples))) {
      stop("Required Stan sampler diagnostics are missing.", call. = FALSE)
    }
    energy <- samples[, "energy__"]
    energy_variance <- stats::var(energy)
    ebfmi <- if (length(energy) > 1L && all(is.finite(energy)) &&
                 is.finite(energy_variance) && energy_variance > 0) {
      mean(diff(energy)^2) / energy_variance
    } else NA_real_
    data.frame(chain = chain,
               divergences = sum(samples[, "divergent__"]),
               treedepth_hits = sum(samples[, "treedepth__"] >= treedepth_limits[chain]),
               ebfmi = ebfmi)
  })
  sampler <- do.call(rbind, rows)
  passed <- all(is.finite(as.matrix(sampler))) &&
    all(sampler$divergences == 0) && all(sampler$treedepth_hits == 0) &&
    all(sampler$ebfmi >= 0.3)
  list(sampler = sampler, passed = isTRUE(passed))
}

diagnostics <- function(fit) {
  .fvc_require_packages()
  draw_array <- rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE)
  parameter_result <- .fvc_parameter_diagnostics(draw_array)
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  limits <- vapply(fit@stan_args, function(arguments) {
    limit <- arguments$control$max_treedepth
    if (is.null(limit)) 10L else as.integer(limit)
  }, integer(1))
  sampler_result <- .fvc_sampler_diagnostics(sampler, limits)
  list(parameters = parameter_result$parameters, sampler = sampler_result$sampler,
       passed = parameter_result$passed && sampler_result$passed)
}
