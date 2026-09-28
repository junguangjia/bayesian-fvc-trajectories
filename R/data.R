# All patient-level data remain private; this module never prints records or IDs.
validate_data <- function(raw) {
  needed <- c("Patient", "Weeks", "FVC", "Age", "Sex", "SmokingStatus")
  if (!is.data.frame(raw) || !all(needed %in% names(raw)))
    stop("Input must contain Patient, Weeks, FVC, Age, Sex and SmokingStatus.")
  if (nrow(raw) < 2L || anyNA(raw[, needed])) stop("Input is empty or contains missing values.")
  if (any(!nzchar(trimws(as.character(raw$Patient))))) stop("Patient keys must be nonempty.")
  for (column in c("Weeks", "FVC", "Age")) {
    if (!is.numeric(raw[[column]]) || any(!is.finite(raw[[column]])))
      stop(paste("Expected finite numeric values for", column))
  }
  if (any(raw$FVC <= 0) || any(raw$Age <= 0)) stop("FVC and Age must be positive.")
  if (!all(raw$Sex %in% c("Female", "Male"))) stop("Unknown Sex category.")
  if (!all(raw$SmokingStatus %in% c("Ex-smoker", "Never smoked", "Currently smokes")))
    stop("Unknown SmokingStatus category.")
  if (anyDuplicated(raw[, needed])) stop("Exact duplicate analysis records require source review.")
  groups <- split(seq_len(nrow(raw)), as.character(raw$Patient))
  for (column in c("Age", "Sex", "SmokingStatus")) {
    if (any(vapply(groups, function(i) length(unique(raw[[column]][i])) != 1L, logical(1))))
      stop(paste("Inconsistent patient-level", column))
  }
  invisible(TRUE)
}

prepare_data <- function(raw, aggregate = TRUE) {
  validate_data(raw)
  raw$Patient <- as.character(raw$Patient)
  raw <- raw[order(raw$Patient, raw$Weeks, raw$FVC), ]
  key <- paste(raw$Patient, format(raw$Weeks, digits = 17), sep = "\r")
  groups <- split(seq_len(nrow(raw)), factor(key, levels = unique(key)))
  if (aggregate) {
    selected <- vapply(groups, `[`, integer(1), 1L)
    d <- raw[selected, , drop = FALSE]
    d$FVC <- vapply(groups, function(i) mean(raw$FVC[i]), numeric(1))
    repeats <- lengths(groups)
  } else {
    d <- raw
    repeats <- rep(1L, nrow(d))
  }
  obs <- data.frame(patient = d$Patient, week = d$Weeks, time = d$Weeks / 52,
    y = d$FVC / 1000, age10 = (d$Age - 65) / 10,
    male = as.integer(d$Sex == "Male"),
    current = as.integer(d$SmokingStatus == "Currently smokes"),
    never = as.integer(d$SmokingStatus == "Never smoked"),
    smoking = d$SmokingStatus, repeat_count = repeats, stringsAsFactors = FALSE)
  obs$id <- match(obs$patient, sort(unique(obs$patient)))
  patients <- obs[!duplicated(obs$patient), c("patient", "id", "age10", "male", "current", "never", "smoking")]
  rownames(obs) <- rownames(patients) <- NULL
  audit <- list(raw_rows = nrow(raw), patients = nrow(patients), patient_weeks = length(groups),
    analyzed_rows = nrow(obs), repeated_week_extra_rows = nrow(raw) - length(groups),
    current_smokers = sum(patients$current), female_current_smokers = sum(patients$current * (1 - patients$male)),
    male_patients = sum(patients$male), female_patients = sum(1 - patients$male),
    ex_smokers = sum(patients$smoking == "Ex-smoker"), never_smokers = sum(patients$never),
    min_week = min(obs$week), max_week = max(obs$week))
  list(obs = obs, patients = patients, audit = audit, aggregate = aggregate)
}

subset_prepared <- function(prepared, patients) {
  out <- prepared
  out$obs <- prepared$obs[prepared$obs$patient %in% patients, , drop = FALSE]
  out$obs$id <- match(out$obs$patient, sort(unique(out$obs$patient)))
  out$patients <- out$obs[!duplicated(out$obs$patient), c("patient", "id", "age10", "male", "current", "never", "smoking")]
  rownames(out$obs) <- rownames(out$patients) <- NULL
  out$audit <- list(patients = nrow(out$patients), analyzed_rows = nrow(out$obs))
  out
}

design_matrix <- function(obs, model) {
  if (!model %in% c("P", "H", "C")) stop("Model must be P, H or C.")
  X <- cbind(alpha = 1, beta_time = obs$time)
  if (model == "C") {
    v <- as.matrix(obs[, c("age10", "male", "current", "never")])
    colnames(v) <- c("beta_age", "beta_male", "beta_current", "beta_never")
    interactions <- v * obs$time
    colnames(interactions) <- c("gamma_age", "gamma_male", "gamma_current", "gamma_never")
    X <- cbind(X, v, interactions)
  }
  X
}

patient_folds <- function(prepared, k = 5L, seed = 4224L) {
  set.seed(seed)
  p <- prepared$patients[order(prepared$patients$patient), ]
  fold <- integer(nrow(p))
  for (group in sort(unique(p$smoking))) {
    i <- which(p$smoking == group)
    fold[i[sample.int(length(i))]] <- rep(seq_len(k), length.out = length(i))
  }
  if (!all(seq_len(k) %in% fold)) stop("Not enough participants to form all folds.")
  data.frame(patient = p$patient, fold = fold, stringsAsFactors = FALSE)
}

simulate_data <- function(n_patients = 50L, seed = 4224L) {
  set.seed(seed)
  age <- round(runif(n_patients, 50, 85))
  sex <- sample(c("Female", "Male"), n_patients, replace = TRUE)
  smoking <- rep(c("Ex-smoker", "Never smoked", "Currently smokes"), length.out = n_patients)
  z <- matrix(rnorm(2 * n_patients), ncol = 2)
  u <- z %*% chol(matrix(c(.6^2, -.2 * .6 * .2, -.2 * .6 * .2, .2^2), 2))
  rows <- lapply(seq_len(n_patients), function(j) {
    weeks <- seq(0, 104, length.out = 8)
    age10 <- (age[j] - 65) / 10
    male <- as.integer(sex[j] == "Male")
    mu <- 3 - .1 * age10 + .4 * male + u[j, 1] +
      (-.25 + .03 * age10 - .04 * male + u[j, 2]) * weeks / 52
    data.frame(Patient = sprintf("SYNTHETIC-%03d", j), Weeks = weeks,
      FVC = (mu + rnorm(8, 0, .12)) * 1000,
      Age = age[j], Sex = sex[j], SmokingStatus = smoking[j])
  })
  do.call(rbind, rows)
}
