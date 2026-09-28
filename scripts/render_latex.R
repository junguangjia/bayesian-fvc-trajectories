#!/usr/bin/env Rscript
# Generate a standalone, aggregate-only LaTeX report. Compilation is separate.
if (!exists("run_dir") || !exists("manifest")) stop("Use scripts/render.R to validate summaries before generating LaTeX.")

read_public <- function(name, columns, required = FALSE) {
  path <- file.path(run_dir, paste0(name, ".csv"))
  if (!file.exists(path)) {
    if (required) stop("Missing public summary: ", name)
    return(as.data.frame(setNames(rep(list(character()), length(columns)), columns)))
  }
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(columns %in% names(x))) stop("Invalid public summary schema: ", name)
  x
}
cohort <- read_public("cohort", c("metric", "value"), TRUE)
estimates <- read_public("estimates", c("model", "variant", "estimand", "mean", "median", "lower95", "upper95", "p_negative"))
diagnostics <- read_public("diagnostics", c("run", "model", "variant", "max_rhat", "min_bulk_ess", "min_tail_ess", "max_slope_mcse", "divergences", "treedepth_hits", "min_ebfmi", "passed"))
cv <- read_public("cv_metrics", c("model", "mean_log_score", "rmse", "mae", "coverage95", "width95"))
cv_diff <- read_public("cv_comparisons", c("model_a", "model_b", "difference", "se"))
check_cols <- c("model", "statistic", "observed", "lower", "median", "upper", "ppp")
checks <- read_public("checks", check_cols)
prior_checks <- read_public("prior_checks", check_cols)
metadata_file <- file.path(run_dir, "fit_metadata.json")
fit_metadata <- if (file.exists(metadata_file)) jsonlite::fromJSON(metadata_file, simplifyVector = FALSE) else list()

scalar <- function(x, fallback = "not recorded") {
  if (is.null(x) || length(x) != 1L || is.na(x)) fallback else as.character(x)
}
# Escape by character so newly inserted LaTeX escape sequences are never escaped again.
tex <- function(x) {
  escapes <- c("\\" = "\\textbackslash{}", "{" = "\\{", "}" = "\\}",
    "#" = "\\#", "$" = "\\$", "%" = "\\%", "&" = "\\&",
    "_" = "\\_", "^" = "\\textasciicircum{}", "~" = "\\textasciitilde{}")
  vapply(as.character(x), function(s) {
    if (is.na(s)) return("not available")
    chars <- strsplit(enc2utf8(s), "", fixed = TRUE)[[1L]]
    replacement <- unname(escapes[chars])
    replacement[is.na(replacement)] <- chars[is.na(replacement)]
    paste0(replacement, collapse = "")
  }, character(1L), USE.NAMES = FALSE)
}
fmt <- function(x, digits = 2L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x)) return("not available")
  formatC(x, format = "f", digits = digits)
}
metric <- function(name) {
  x <- cohort$value[cohort$metric == name]
  if (length(x) != 1L) return(NA_real_)
  suppressWarnings(as.numeric(x))
}
is_true <- function(x) !is.na(x) & tolower(as.character(x)) %in% c("true", "1")
diagnostics_ok <- nrow(diagnostics) > 0L && all(is_true(diagnostics$passed))
is_research <- identical(scalar(manifest$mode), "research")
scientific_ok <- is_research && identical(scalar(manifest$status), "verified") &&
  diagnostics_ok && nrow(diagnostics) >= 24L
primary <- estimates[estimates$variant == "primary", , drop = FALSE]
main_slopes <- primary[primary$estimand == "standardized_slope", , drop = FALSE]
h_main <- main_slopes[main_slopes$model == "H", , drop = FALSE]
p_main <- main_slopes[main_slopes$model == "P", , drop = FALSE]
ch <- cv_diff[cv_diff$model_a == "C" & cv_diff$model_b == "H", , drop = FALSE]
serial <- checks[checks$model %in% c("H", "C") & checks$statistic == "residual_lag1_correlation", , drop = FALSE]
flagged <- serial[is.finite(serial$ppp) & (serial$ppp < .025 | serial$ppp > .975), , drop = FALSE]
association_labels <- c(age_main_ml_per_decade = "Age, per decade", male_main_ml = "Male vs female",
  current_main_ml = "Current vs former smoker", never_main_ml = "Never vs former smoker",
  age_time_ml_week_per_decade = "Age, per decade", male_time_ml_week = "Male vs female",
  current_time_ml_week = "Current vs former smoker", never_time_ml_week = "Never vs former smoker")
profile_labels <- c(reference_slope = "Female, former smoker, age 65",
  male_ex_smoker_age65_slope = "Male, former smoker, age 65",
  female_never_smoker_age65_slope = "Female, never smoker, age 65",
  female_current_smoker_age65_slope = "Female, current smoker, age 65",
  female_ex_smoker_age75_slope = "Female, former smoker, age 75")
main_labels <- c(standardized_slope = "Overall slope (mL/week)", tau_intercept = "Intercept SD (mL)",
  tau_slope = "Slope SD (mL/week)", rho = "Intercept--slope correlation")
check_labels <- c(negative_fraction = "Negative FVC fraction", median_fvc_ml = "FVC median (mL)",
  q01_fvc_ml = "FVC 1st percentile (mL)", q99_fvc_ml = "FVC 99th percentile (mL)",
  median_mL = "FVC median (mL)", q05_mL = "FVC 5th percentile (mL)",
  q95_mL = "FVC 95th percentile (mL)", sd_mL = "FVC SD (mL)",
  patient_mean_sd_mL = "SD of patient means (mL)",
  patient_slope_mean_mL_per_year = "Mean patient slope (mL/year)",
  patient_slope_sd_mL_per_year = "SD of patient slopes (mL/year)",
  residual_lag1_correlation = "Adjacent residual correlation",
  residual_time_slope_mL_per_year = "Residual time slope (mL/year)",
  absolute_residual_time_slope_mL_per_year = "Absolute-residual time slope (mL/year)")
label <- function(x, mapping) {
  result <- unname(mapping[as.character(x)])
  result[is.na(result)] <- as.character(x[is.na(result)])
  result
}
variant_label <- function(x) gsub("_", "-", sub("cv_fold_", "CV", x, fixed = TRUE), fixed = TRUE)

lines <- character()
emit <- function(...) lines <<- c(lines, paste0(...))
paragraph <- function(...) { emit(...); emit("") }
row <- function(cells) emit(paste(cells, collapse = " & "), " \\\\")
table_start <- function(spec, headers, caption) {
  spec <- gsub("p{", ">{\\raggedright\\arraybackslash}p{", spec, fixed = TRUE)
  emit("\\begingroup\\small", "\\setlength{\\tabcolsep}{4pt}",
       paste0("\\begin{longtable}{", spec, "}"),
       paste0("\\caption{", tex(caption), "}\\\\"), "\\toprule")
  row(tex(headers)); emit("\\midrule\\endfirsthead", "\\toprule")
  row(tex(headers)); emit("\\midrule\\endhead", "\\bottomrule\\endfoot")
}
table_end <- function() emit("\\end{longtable}", "\\endgroup", "")
interval <- function(lo, hi, digits = 2L) paste0("[", fmt(lo, digits), ", ", fmt(hi, digits), "]")
empty_note <- function() paragraph("This component is unavailable in the selected run.")
forest <- function(x, labels, caption) {
  valid <- is.finite(x$mean) & is.finite(x$lower95) & is.finite(x$upper95)
  x <- x[valid, , drop = FALSE]; labels <- labels[valid]
  if (!nrow(x)) return(invisible(NULL))
  n <- nrow(x); limits <- range(c(0, x$lower95, x$upper95))
  pad <- max(diff(limits) * .08, .05)
  ticks <- pretty(limits + c(-pad, pad), n = 5L)
  limits <- range(ticks)
  left <- 110; plot_width <- 325; axis_y <- 32; height <- 45 + 23 * n
  xpos <- function(value) left + plot_width * (value - limits[1]) / diff(limits)
  point <- function(value) formatC(value, format = "f", digits = 3L)
  emit("\\begin{figure}[htbp]\\centering\\setlength{\\unitlength}{1pt}")
  emit("\\begin{picture}(455,", height, ")")
  emit("\\put(", left, ",", axis_y, "){\\line(1,0){", plot_width, "}}")
  emit("\\multiput(", point(xpos(0)), ",", axis_y, ")(0,6){",
       floor((height - axis_y - 8) / 6), "}{\\line(0,1){2}}")
  tick_digits <- if (diff(limits) < 2) 2L else if (diff(limits) < 10) 1L else 0L
  for (tick in ticks) {
    emit("\\put(", point(xpos(tick)), ",", axis_y - 2, "){\\line(0,1){4}}")
    emit("\\put(", point(xpos(tick)), ",20){\\makebox(0,0){\\scriptsize $", fmt(tick, tick_digits), "$}}")
  }
  emit("\\put(", left + plot_width / 2, ",5){\\makebox(0,0){\\small Slope (mL/week)}}")
  for (i in seq_len(n)) {
    y <- 54 + 23 * (n - i)
    lo <- xpos(x$lower95[i]); hi <- xpos(x$upper95[i])
    emit("\\put(", left - 8, ",", y, "){\\makebox(0,0)[r]{\\small ", tex(labels[i]), "}}")
    emit("\\put(", point(lo), ",", y, "){\\line(1,0){", point(hi - lo), "}}")
    emit("\\put(", point(lo), ",", y - 3, "){\\line(0,1){6}}")
    emit("\\put(", point(hi), ",", y - 3, "){\\line(0,1){6}}")
    emit("\\put(", point(xpos(x$mean[i])), ",", y, "){\\circle*{4}}")
  }
  emit("\\end{picture}\\caption{", tex(caption), "}\\end{figure}")
  emit("")
}
check_table <- function(x, caption) {
  if (!nrow(x)) return(empty_note())
  table_start("lp{0.32\\textwidth}rrr", c("Model", "Statistic", "Observed", "95% replicated", "Tail fraction"), caption)
  for (i in seq_len(nrow(x))) {
    z <- x[i, ]; digits <- if (grepl("correlation|fraction", z$statistic)) 3L else 1L
    row(c(tex(z$model), tex(label(z$statistic, check_labels)), fmt(z$observed, digits),
          interval(z$lower, z$upper, digits), fmt(z$ppp, 3L)))
  }
  table_end()
}

emit("\\documentclass[11pt,letterpaper]{article}",
  "\\usepackage[margin=1in]{geometry}", "\\usepackage[T1]{fontenc}",
  "\\usepackage{lmodern}", "\\usepackage{amsmath,amssymb,booktabs,longtable,array}",
  "\\usepackage[hidelinks]{hyperref}", "\\hypersetup{pdftitle={Bayesian Modeling of Pulmonary Fibrosis Trajectories},pdfauthor={Junguang Jia}}", "\\setlength{\\parskip}{0.35em}",
  "\\setlength{\\parindent}{0pt}", "\\setlength{\\emergencystretch}{2em}",
  "\\title{Bayesian Modeling of Pulmonary Fibrosis Trajectories}",
  "\\author{Junguang Jia}", "\\date{}", "\\begin{document}", "\\maketitle")
paragraph("\\textbf{Run:} ", tex(run_id), "\\quad\\textbf{Mode:} ", tex(scalar(manifest$mode)),
          "\\quad\\textbf{Status:} ", tex(scalar(manifest$status)))
if (scientific_ok) {
  paragraph("\\textbf{Computational status.} The required research fits meet the recorded sampling gates. This establishes computational acceptance, not substantive model adequacy; predictive discrepancies and observational limitations remain relevant.")
} else if (!is_research) {
  paragraph("\\textbf{Synthetic demonstration / software check.} All numerical summaries concern generated data. They are not findings about pulmonary fibrosis. Reduced smoke-test sampling does not establish convergence.")
} else {
  paragraph("\\textbf{Unverified research output.} Required fits or computational checks are incomplete. Available summaries are retained for diagnosis and must not be presented as validated scientific conclusions.")
}
emit("\\begin{abstract}")
paragraph("This exploratory longitudinal study estimates forced vital capacity (FVC) change and between-patient variation, examines covariate associations, and assesses prediction for patients without observed FVC histories. Pooled, hierarchical, and covariate-adjusted Gaussian models are compared using patient-level validation. The selected input contains ",
  fmt(metric("patients"), 0), " patients, ", fmt(metric("raw_rows"), 0), " measurements, and ",
  fmt(metric("patient_weeks"), 0), " distinct patient-weeks.")
if (scientific_ok && nrow(h_main) == 1L) {
  paragraph("The primary hierarchical posterior mean slope is ", fmt(h_main$mean),
    " mL/week (95\\% credible interval ", interval(h_main$lower95, h_main$upper95),
    "). This is a model-based observational estimate, not a causal effect. ",
    if (nrow(flagged)) "Residual predictive checks flag remaining serial dependence, limiting model adequacy and interval interpretation." else
      "Sensitivity and predictive checks accompany the main estimate.")
} else paragraph("No validated scientific effect estimate is asserted for this run.")
if (scientific_ok && nrow(ch) == 1L) paragraph("In internal patient-level validation, the adjusted model's mean joint log score exceeds H's by ", fmt(ch$difference, 3), " (descriptive SE ", fmt(ch$se, 3), "). This comparison concerns previously unseen patients from the same data source and does not establish external validity.")
emit("\\end{abstract}", "\\section{Introduction}")
paragraph("Longitudinal FVC provides information about overall change and heterogeneity in patient trajectories. The primary question concerns a population time slope. Adjusted associations and prediction for a previously unseen patient are distinct secondary questions. This analysis is exploratory and does not evaluate treatment effects or a clinical decision rule.")
paragraph("FVC is a spirometric measure of lung function. The OSIC competition describes a prediction task using baseline imaging and clinical information~\\cite{osic}. The present study uses only tabular covariates and asks about a patient with no observed FVC history; it is not a reproduction of that competition task or a comparison with leaderboard scores.")
emit("\\clearpage", "\\section{Data}")
table_start("lr", c("Aggregate characteristic", "Value"), "Cohort summaries calculated from this run's input.")
for (i in seq_len(nrow(cohort))) row(c(tex(gsub("_", " ", cohort$metric[i], fixed = TRUE)), tex(cohort$value[i])))
table_end()
paragraph("The main analysis averages FVC for repeated records from the same patient and week; an H/C sensitivity analysis retains every original row. Invalid inputs, missing model fields, inconsistent patient-level covariates, and exact duplicate analysis records require source review. The Percent field is excluded. Original week zero is retained and is not interpreted as diagnosis.")
if (is.finite(metric("current_smokers"))) paragraph("There are ", fmt(metric("current_smokers"), 0),
  " current smokers, including ", fmt(metric("female_current_smokers"), 0),
  " female current smokers. Associations involving these small subgroups are exploratory.")
paragraph("The research input schema resembles the tabular component of the \\href{https://www.kaggle.com/competitions/osic-pulmonary-fibrosis-progression/data}{OSIC Pulmonary Fibrosis Progression competition}; exact provenance of the local file has not been independently verified. Obtain authorized data directly from its provider and review the \\href{https://www.kaggle.com/competitions/osic-pulmonary-fibrosis-progression/rules}{applicable rules}. No real patient records, identifiers, individual predictions, or observed trajectories are distributed. Synthetic demonstrations are independently generated.")
emit("\\section{Methods}", "\\subsection{Models and estimands}")
paragraph("For patient $i$ at visit $j$, set $y_{ij}=\\mathrm{FVC}_{ij}/1000$ (L), $t_{ij}=\\mathrm{Weeks}_{ij}/52$ (years), and $a_i=(\\mathrm{Age}_i-65)/10$. Let $x_i$ contain centered age, male, current-smoker, and never-smoker indicators. The reference profile is female, former smoker, age 65.")
emit("\\begin{align*}",
  "\\mathrm{P}:\\quad y_{ij}&=\\alpha+\\beta_t t_{ij}+\\epsilon_{ij},\\\\",
  "\\mathrm{H}:\\quad y_{ij}&=\\alpha+\\beta_t t_{ij}+b_{0i}+b_{1i}t_{ij}+\\epsilon_{ij},\\\\",
  "\\mathrm{C}:\\quad y_{ij}&=\\alpha+\\beta_t t_{ij}+x_i^\\top\\gamma+t_{ij}x_i^\\top\\delta",
  "+b_{0i}+b_{1i}t_{ij}+\\epsilon_{ij},\\\\",
  "\\epsilon_{ij}&\\sim N(0,\\sigma^2),\\qquad b_i\\sim N_2(0,D\\Omega D),",
  "\\end{align*}")
paragraph("Here $D=\\operatorname{diag}(\\tau_0,\\tau_1)$. H and C use the noncentered representation $b_i=DL_\\Omega z_i$, $z_i\\sim N_2(0,I)$, with no sample sum-to-zero constraint. Conditional residuals are independent. The \\href{https://mc-stan.org/docs/stan-users-guide/regression.html}{Stan User's Guide} describes this hierarchical construction.")
emit("\\[s_H=\\frac{1000}{52}\\beta_t,\\qquad{}",
  "s_C=\\frac{1000}{52}\\left(\\beta_t+\\frac1N\\sum_{i=1}^N x_i^\\top\\delta\\right).\\]")
paragraph("The C overall slope gives each patient's covariates equal weight. Its reference-profile slope is $1000\\beta_t/52$, a separate estimand. All principal slopes and slope SDs are reported in mL/week; intercept SDs and level contrasts are in mL.")
emit("\\subsection{Priors and sampling}")
table_start("ll", c("Parameter", "Main prior on transformed scale"), "Analyst-specified priors; locations and scales are not estimated from FVC outcomes.")
prior_rows <- list(c("Intercept", "$N(3,1.5^2)$"), c("Time slope", "$N(0,0.5^2)$"),
  c("Age main effect", "$N(0,0.5^2)$"), c("Sex or smoking main effect", "$N(0,1)$"),
  c("Time by age", "$N(0,0.25^2)$"), c("Time by sex or smoking", "$N(0,0.5^2)$"),
  c("Random-intercept SD", "$\\operatorname{HalfNormal}(0,1)$"),
  c("Random-slope SD; residual SD", "$\\operatorname{HalfNormal}(0,0.5)$"),
  c("Random-effect correlation", "$\\Omega\\sim\\operatorname{LKJ}(2)$"))
for (z in prior_rows) row(c(tex(z[1]), z[2])); table_end()
paragraph("HalfNormal arguments give location and scale. Prior predictive simulation precedes fitting and conditions on the design, not FVC. H/C sensitivity fits multiply every normal and half-normal scale by 0.5 and 2, retaining locations and LKJ(2). Gaussian priors may imply negative FVC; this is exposed in the prior checks rather than silently truncated.")
paragraph("The research workflow contains 24 fits: three main fits, 15 validation fits, four prior-scale sensitivity fits, and two all-row sensitivity fits. Defaults are four chains, 6,000 iterations including 3,000 warmup, adapt-delta 0.95, maximum tree depth 12, and base seed 4224. A failed non-smoke diagnostic gate permits one retry at 10,000 iterations, 5,000 warmup, and adapt-delta 0.99; maximum depth becomes 14 only after tree-depth hits. If that retry fails, a targeted third attempt raises maximum depth to 14 for depth hits and extends to 20,000 iterations with 5,000 warmup for R-hat or ESS failure, retaining adapt-delta 0.99. A further failure stops the run as unverified. These computational extensions leave data, priors, and model structure unchanged and do not use predictive performance. Original attempts are preserved; actual accepted settings appear below.")
emit("\\subsection{New-patient validation}")
paragraph("Five folds are assigned by patient, stratified by smoking status using seed 4224. P/H/C share the allocation. Held-out FVC is used only for scoring. Covariates and visit times are available, but no patient's held-out measurements estimate their random effects. Conditional on one posterior draw, the joint response is")
emit("\\[y_i\\mid\\theta,\\mathcal D_{-k}\\sim N_{n_i}\\!\\left(X_i\\beta,",
  "Z_iD\\Omega DZ_i^\\top+\\sigma^2 I\\right),\\qquad Z_i=[\\mathbf1,t_i].\\]")
paragraph("P omits the random-effect covariance. For each patient, the joint predictive density is averaged over all retained posterior draws before taking its logarithm. The reported log score is on the mL scale, using the adjustment $-n_i\\log(1000)$ to the litre-scale log density. It is not a sum of separately integrated visit scores. The primary aggregate gives each patient equal weight.")
paragraph("RMSE is the square root of the patient-average within-patient MSE. MAE, interval coverage, and width likewise first average over visits within each patient. Up to 1,000 posterior predictive draws preserve shared patient effects across visits. Intervals are marginal 95\\% visit intervals, not simultaneous bands. Paired score-difference standard errors are descriptive and do not fully reflect dependence across folds. The prediction target follows the group-level distinction in the \\href{https://mc-stan.org/loo/articles/online-only/faq.html}{loo cross-validation FAQ}; this is internal, not external, validation.")

emit("\\section{Results}", "\\subsection{Population change and heterogeneity}")
z <- primary[primary$estimand %in% names(main_labels), , drop = FALSE]
if (nrow(z)) {
  table_start("llrrr", c("Model", "Estimand", "Mean", "95% credible interval", "Pr(negative)"), "Primary posterior summaries. Correlations are unitless.")
  for (i in seq_len(nrow(z))) row(c(tex(z$model[i]), tex(label(z$estimand[i], main_labels)),
    fmt(z$mean[i]), interval(z$lower95[i], z$upper95[i]), fmt(z$p_negative[i], 3L)))
  table_end()
} else empty_note()
forest(main_slopes, paste0(main_slopes$model, ": overall slope"), "Posterior means and 95% credible intervals for standardized overall slopes. These remain subject to the run's stated verification status.")
paragraph("Pr(negative) is estimated from retained posterior draws, not a frequentist p-value. Values rounded to 1.000 do not imply literal certainty. For positive scale parameters it is not an inferential test. The random-effect correlation refers to the original week-zero intercept and can depend on that time origin.")
emit("\\clearpage", "\\subsection{Exploratory covariate associations}")
z <- primary[primary$model == "C" & primary$estimand %in% names(association_labels), , drop = FALSE]
if (nrow(z)) {
  table_start("p{0.32\\textwidth}lrr", c("Contrast", "Quantity", "Mean", "95% credible interval"), "Adjusted C contrasts. Age contrasts are per decade; level contrasts apply at original week zero.")
  for (i in seq_len(nrow(z))) row(c(tex(label(z$estimand[i], association_labels)),
    if (grepl("_main_", z$estimand[i])) "Level (mL)" else "Slope (mL/week)",
    fmt(z$mean[i]), interval(z$lower95[i], z$upper95[i])))
  table_end()
} else empty_note()
z <- primary[primary$model == "C" & primary$estimand %in% names(profile_labels), , drop = FALSE]
if (nrow(z)) {
  table_start("p{0.48\\textwidth}rr", c("Profile", "Mean", "95% credible interval"), "Selected C profile slopes (mL/week), separate from the cohort-standardized slope.")
  for (i in seq_len(nrow(z))) row(c(tex(label(z$estimand[i], profile_labels)), fmt(z$mean[i]), interval(z$lower95[i], z$upper95[i])))
  table_end()
}
paragraph("These contrasts hold other modeled covariates fixed. They describe associations, not causal effects. Small subgroups warrant particular caution even when Monte Carlo error is small.")
emit("\\subsection{Sensitivity analyses}")
z <- estimates[estimates$model %in% c("H", "C") & estimates$estimand == "standardized_slope", , drop = FALSE]
if (any(z$variant != "primary")) {
  table_start("llrr", c("Model", "Variant", "Mean", "95% credible interval"), "Overall slopes under prior-scale and duplicate-handling variants (mL/week).")
  for (i in seq_len(nrow(z))) row(c(tex(z$model[i]), tex(variant_label(z$variant[i])), fmt(z$mean[i]), interval(z$lower95[i], z$upper95[i])))
  table_end()
  forest(z, paste(z$model, variant_label(z$variant), sep = ": "), "Sensitivity of H and C overall slopes: posterior means and 95% credible intervals.")
} else empty_note()
paragraph("Prior-half and prior-double change only the prespecified prior scales; all-rows retains every original measurement. Agreement across these variants does not establish robustness to other assumptions.")

emit("\\section{Validation}", "\\subsection{Sampling diagnostics}")
paragraph("Computational gates require rank-normalized R-hat below 1.01, bulk/tail ESS at least 400, no post-warmup divergences, and each chain's E-BFMI at least 0.3 for monitored parameters and major derived estimands. These follow \\href{https://mc-stan.org/learn-stan/diagnostics-warnings.html}{Stan's diagnostic guidance} and rank-based convergence diagnostics~\\cite{rhat}. This project additionally requires zero tree-depth hits. Slope MCSE is reported in mL/week. Passing these checks does not establish model adequacy.")
if (nrow(diagnostics)) {
  table_start("lrrrrrrrr", c("Fit", "R-hat", "Bulk", "Tail", "Div.", "Depth", "BFMI", "MCSE", "Pass"), "Diagnostic extrema for each accepted or most recently attempted fit.")
  for (i in seq_len(nrow(diagnostics))) {
    z <- diagnostics[i, ]
    row(c(tex(paste(z$model, variant_label(z$variant))), fmt(z$max_rhat, 4),
      fmt(z$min_bulk_ess, 0), fmt(z$min_tail_ess, 0), fmt(z$divergences, 0),
      fmt(z$treedepth_hits, 0), fmt(z$min_ebfmi, 2), fmt(z$max_slope_mcse, 3),
      if (is_true(z$passed)) "yes" else "no"))
  }
  table_end()
} else empty_note()
emit("\\subsection{Patient-level predictive performance}")
if (nrow(cv)) {
  table_start("lrrrrr", c("Model", "Mean log score", "RMSE", "MAE", "Coverage", "Width"), "Patient-equal validation: mL-scale joint log density; errors and width in mL; coverage as a proportion.")
  for (i in seq_len(nrow(cv))) {
    z <- cv[i, ]; row(c(tex(z$model), fmt(z$mean_log_score, 3), fmt(z$rmse, 1),
      fmt(z$mae, 1), fmt(z$coverage95, 3), fmt(z$width95, 1)))
  }
  table_end()
} else empty_note()
if (nrow(cv_diff)) {
  table_start("lrr", c("Comparison", "Mean paired difference", "Approximate SE"), "Joint-log-score contrasts, model A minus B. Positive values favor A on this internal criterion.")
  for (i in seq_len(nrow(cv_diff))) row(c(tex(paste(cv_diff$model_a[i], "minus", cv_diff$model_b[i])), fmt(cv_diff$difference[i], 3), fmt(cv_diff$se[i], 3)))
  table_end()
}
paragraph("Joint scores depend on visit number and timing; models share the same held-out records. Coverage must be interpreted with interval width. Neither a favorable score nor nominal empirical coverage establishes transportability.")
emit("\\subsection{Prior and posterior predictive checks}")
check_table(prior_checks, "Prior checks: central 95% simulated statistic intervals and upper-tail fractions.")
check_table(checks, "Conditional posterior checks. Slopes in this table use mL/year, not mL/week.")
paragraph("Posterior checks condition on observed patients and their posterior random effects. Patient OLS slopes require at least two distinct times. Adjacent residual pairs are pooled within patients after ordering by time; equal spacing is not assumed. Signed and absolute residual time slopes address different discrepancies. Distributional observed statistics are fixed; observed residual entries are medians of draw-specific discrepancies. Tail fractions compare replicated and observed discrepancies within the same posterior draw before averaging. Values below 0.025 or above 0.975 flag discrepancies, not calibrated hypothesis-test significance.")

emit("\\section{Discussion}")
if (scientific_ok) paragraph("The main estimates have passed the recorded computational gates. Substantive interpretation remains conditional on the likelihood, priors, observation process, and predictive-check findings.") else paragraph("This run does not support a validated research conclusion. The available numerical output describes a synthetic demonstration or an incomplete research analysis.")
if (scientific_ok && nrow(p_main) == 1L && p_main$lower95 < 0 && p_main$upper95 > 0) paragraph("The pooled benchmark's 95\\% slope interval includes zero. P combines within- and between-patient information without patient effects, whereas H is the prespecified primary model for longitudinal change. The decline conclusion therefore should not be described as independently established by all three models.")
if (nrow(serial)) for (i in seq_len(nrow(serial))) {
  z <- serial[i, ]; flag <- is.finite(z$ppp) && (z$ppp < .025 || z$ppp > .975)
  paragraph("For model ", tex(z$model), ", median observed adjacent-visit residual correlation is ",
    fmt(z$observed, 3), ", versus replicated 95\\% interval ", interval(z$lower, z$upper, 3),
    " (simulated upper-tail fraction ", fmt(z$ppp, 3), "). ",
    if (flag) "This flags remaining serial dependence." else "This check is not flagged at the chosen thresholds; absence of a flag does not establish independence.")
}
if (nrow(flagged)) paragraph("Remaining serial dependence makes conditional residual independence questionable. Population credible intervals and new-patient predictive intervals retain that assumption; their calibration and predictive adequacy may be affected. Passing MCMC diagnostics does not resolve this limitation. The prespecified models are retained so the discrepancy remains visible.")
ch <- cv_diff[cv_diff$model_a == "C" & cv_diff$model_b == "H", , drop = FALSE]
if (scientific_ok && nrow(ch) == 1L) paragraph("The paired mean joint-log-score difference C minus H is ",
  fmt(ch$difference, 3), " (approximate SE ", fmt(ch$se, 3),
  "). This is an internal comparison with a descriptive patient-level uncertainty approximation, not evidence of general model superiority.")
if (scientific_ok && nrow(cv)) for (model in c("H", "C")) {
  z <- cv[cv$model == model, , drop = FALSE]
  if (nrow(z) == 1L) paragraph("Model ", tex(model), " has patient-equal empirical coverage ",
    fmt(100 * z$coverage95, 1), "\\% for nominal 95\\% intervals, with mean width ", fmt(z$width95, 1),
    " mL. Coverage and width must be considered together.")
}
h_prior <- estimates[estimates$model == "H" & estimates$estimand == "standardized_slope" & estimates$variant %in% c("primary", "prior_half", "prior_double"), , drop = FALSE]
if (scientific_ok && all(c("primary", "prior_half", "prior_double") %in% h_prior$variant) && all(is.finite(h_prior$mean))) {
  r <- range(h_prior$mean); base <- h_prior$mean[h_prior$variant == "primary"]
  paragraph("Across H's main, half-scale, and double-scale priors, posterior mean slopes range from ",
    fmt(r[1]), " to ", fmt(r[2]), " mL/week. ",
    if (length(base) == 1L) paste0("The largest absolute change from the main estimate is ", fmt(max(abs(h_prior$mean - base))), " mL/week. ") else "",
    "This addresses those specific scale choices, not other prior families or residual dependence.")
}
paragraph("Linear trajectories, Gaussian constant-variance errors, informative visit timing or dropout, cohort selection, and unmeasured factors limit inference. Averaging same-week outcomes changes their weighting and noise structure. Covariate adjustment does not make associations causal. Small subgroups can remain prior sensitive. Population-mean credible intervals and individual predictive intervals answer different questions.")
paragraph("Internal prediction at known covariates and visit times does not establish external validity, future-visit performance for a known patient, treatment benefit, or clinical utility. Independent data with documented provenance and models addressing informative observation processes would be needed for those questions.")

emit("\\clearpage", "\\section{Reproducibility}")
table_start("ll", c("Field", "Recorded value"), "Run provenance without participant-level records or local paths.")
for (field in c("run_id", "mode", "created_utc", "status", "seed", "run_count")) row(c(tex(gsub("_", " ", field, fixed = TRUE)), tex(scalar(manifest[[field]]))))
table_end()
paragraph("Input SHA-256: \\nolinkurl{", tex(scalar(manifest$data_sha256)), "}.")
if (!is.null(manifest$software)) {
  sw <- unlist(manifest$software)
  table_start("ll", c("Component", "Version"), "Software recorded by the run.")
  for (i in seq_along(sw)) row(c(tex(names(sw)[i]), tex(sw[i]))); table_end()
}
if (length(fit_metadata)) {
  table_start("lrrrrrrr", c("Fit", "Chains", "Iter.", "Warmup", "Delta", "Depth", "Seconds", "Attempts"), "Actual fit settings and elapsed times. Attempts include retained initial failures.")
  for (entry in fit_metadata) {
    s <- entry$settings
    row(c(tex(paste(entry$model, variant_label(entry$variant))),
      tex(scalar(s$chains)), tex(scalar(s$iter)), tex(scalar(s$warmup)),
      tex(scalar(s$adapt_delta)), tex(scalar(s$max_treedepth)),
      fmt(entry$elapsed_seconds, 1), as.character(length(entry$attempts))))
  }
  table_end()
}
paragraph("The report reads only the selected run's public summaries and manifest. Cache identity includes data, preprocessing, model, priors, seeds, sampling settings, and software. Stale or incompatible caches cannot silently supply results. Real data, fitted objects, and patient-level validation records remain private.")
emit(c("\\begin{verbatim}", "Rscript scripts/setup.R", "Rscript tests/run_tests.R",
  "Rscript scripts/run.R --mode demo --smoke", "Rscript scripts/render.R --run demo", "",
  "Rscript scripts/run.R --mode research --data <authorized.csv>",
  "Rscript scripts/render.R --run <run-id>", "\\end{verbatim}"))
paragraph("Code and original documentation: \\href{https://github.com/junguangjia/bayesian-fvc-trajectories}{github.com/junguangjia/bayesian-fvc-trajectories}, under the MIT License. That license does not grant third-party data rights. Cite the data provider separately.")
emit("\\begin{thebibliography}{9}\\raggedright",
  "\\bibitem{osic} Open Source Imaging Consortium. OSIC Pulmonary Fibrosis Progression: Dataset Description. Kaggle, 2020. \\url{https://www.kaggle.com/competitions/osic-pulmonary-fibrosis-progression/data}.",
  "\\bibitem{rhat} A. Vehtari, A. Gelman, D. Simpson, B. Carpenter, and P.-C. Buerkner. Rank-normalization, folding, and localization: An improved $\\widehat R$ for assessing convergence of MCMC. Bayesian Analysis, 2021. \\url{https://arxiv.org/abs/1903.08008}.",
  "\\bibitem{loo} Stan Development Team. Frequently asked questions about cross-validation. \\url{https://mc-stan.org/loo/articles/online-only/faq.html}.",
  "\\bibitem{stan} Stan Development Team. Stan User's Guide: Regression models. \\url{https://mc-stan.org/docs/stan-users-guide/regression.html}.",
  "\\end{thebibliography}", "\\end{document}")
writeLines(lines, file.path(run_dir, "report.tex"), useBytes = TRUE)
cat("Standalone LaTeX written to", file.path("results", run_id, "report.tex"), "\n")
