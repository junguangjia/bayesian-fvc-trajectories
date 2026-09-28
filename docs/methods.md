# Statistical methods

## Questions and estimands

The study estimates longitudinal FVC change, quantifies variation between patients, and explores associations with age, sex, and smoking status. It separately evaluates prediction for a new patient whose covariates and measurement times are known but whose FVC has not been observed.

The primary inferential model is H. Its population time coefficient is the mean slope in the model's patient distribution. Model C's overall slope averages its fixed-effect slope over the observed covariate distribution, assigning each patient equal weight. Neither estimand is a measurement-weighted average of fitted patient-specific slopes. The reference-profile slope and other profile slopes are separate estimands.

## Likelihood and hierarchical structure

Let `y_ij = FVC_ij / 1000` denote FVC in litres, `t_ij = Weeks_ij / 52` time in years, and `a_i = (Age_i - 65) / 10` age in decades from 65. Within each model, errors are independent Gaussian variables conditional on its parameters and random effects.

**P:** `y_ij ~ Normal(alpha + beta_t * t_ij, sigma)`.

**H:** `y_ij ~ Normal(alpha + beta_t * t_ij + b_0i + b_1i * t_ij, sigma)`.

**C:** the H mean additionally contains age, sex, and smoking-status main effects and their interactions with time. The reference profile is female, former smoker (`Ex-smoker`), and age 65. Its intercept is at original week zero. These are coding references, not a population average.

For H and C, `(b_0i, b_1i)` has a zero-mean bivariate normal population distribution with covariance `diag(tau) * Omega * diag(tau)`. The Stan implementation uses standard-normal latent variables and a Cholesky factor to represent this distribution without imposing a sum-to-zero constraint across the observed patients. This follows the multivariate hierarchical-regression construction in the [Stan User's Guide](https://mc-stan.org/docs/stan-users-guide/regression.html#multivariate-priors-for-hierarchical-models).

For C, write the fixed slope for patient covariates `x_i` as `s_i = beta_t + x_i' * beta_interaction`. Each posterior draw of the standardized overall slope is `mean_i(s_i)`. Multiply a slope in litres/year by `1000 / 52` to report mL/week. This conversion also applies to the slope's interval bounds and Monte Carlo standard error.

## Prespecified priors

All priors are analyst specified on the transformed scale, independently of the supplied FVC mean and standard deviation.

| Quantity | Main prior |
|:--|:--|
| Intercept (L) | Normal(3, 1.5) |
| Time slope (L/year) | Normal(0, 0.5) |
| Age main effect (L/decade) | Normal(0, 0.5) |
| Sex and smoking main effects (L) | Normal(0, 1) |
| Time by age (L/year/decade) | Normal(0, 0.25) |
| Time by sex or smoking (L/year) | Normal(0, 0.5) |
| Random-intercept standard deviation (L) | HalfNormal(0, 1) |
| Random-slope standard deviation (L/year) | HalfNormal(0, 0.5) |
| Residual standard deviation (L) | HalfNormal(0, 0.5) |
| Random-effect correlation matrix | LKJ(2) |

Prior predictive simulation occurs before posterior fitting. It conditions on the available design values, but not on FVC outcomes. The report compares broad FVC distribution summaries under these priors with the observed summaries. The Gaussian model and broad intercept prior can generate physiologically impossible negative FVC; their frequency and extent are a model-checking issue, not an implicit truncation of the sampling model.

Sensitivity fits multiply every normal or half-normal prior scale in H and C by 0.5 and by 2, retaining prior locations and the LKJ(2) correlation prior. Hyperparameters are not optimized against held-out prediction results.

## New-patient validation

Five folds are assigned at the patient level, stratified by smoking status with seed 4224. The same allocation is used for all three models. All records from a held-out patient remain out of the training data. Covariates and requested visit times are permitted prediction inputs; held-out FVC values are used only for scoring.

For a held-out patient with design matrix `X_i`, time matrix `Z_i = [1, t_i]`, and one posterior draw of population parameters, the entire response vector follows:

`y_i | theta, training ~ MVN(X_i * beta, Z_i * Sigma_b * Z_i' + sigma^2 * I)`.

For P, omit the random-effect covariance term. The primary score integrates this joint density over posterior draws using log-mean-exp. It is a joint patient score, not the sum of separately integrated visit scores. This preserves dependence induced by common patient effects and uncertainty in shared parameters.

The primary aggregate is the mean joint log predictive density across patients, evaluated using all retained posterior draws. Scores are reported on the mL outcome scale: the joint log density computed in litres is adjusted by subtracting `n_i * log(1000)` for a patient with `n_i` visits. Changing units shifts all models' scores for a given patient by the same constant. Pairwise comparisons use matched patient score differences. Their reported standard error is the sample standard deviation of those differences divided by the square root of the patient count; it is descriptive and does not account fully for dependence across folds or uncertainty in the fold assignment.

Auxiliary errors and coverage weight patients equally: average each patient's visit-level squared error, absolute error, interval coverage, or interval width before averaging across patients. RMSE is the square root of the resulting mean squared error. Predictive intervals come from joint posterior predictive simulations sharing one random-effect vector within a new patient's visits. The reported intervals are 95% marginal intervals at each visit, not simultaneous trajectory bands. Predictive simulation uses at most 1,000 retained posterior draws.

The target is a new patient from the same data source, rather than a future visit for an already observed patient. The [loo cross-validation FAQ](https://mc-stan.org/loo/articles/online-only/faq.html) explains why the held-out unit must match the prediction target. This design supports internal validation only.

## Computation, checks, and sensitivity

The full research workflow includes 24 fits: P/H/C on the main data, P/H/C in each of five folds, four prior-sensitivity H/C fits, and H/C fits retaining all duplicate measurements. Default sampling is four chains, 6,000 iterations including 3,000 warmup, `adapt_delta = 0.95`, and `max_treedepth = 12`. Actual settings and elapsed time are recorded for each fit.

Sampling acceptance requires rank-normalized R-hat below 1.01 and bulk/tail effective sample sizes of at least 400 for monitored parameters and main derived estimands, no post-warmup divergences, and E-BFMI of at least 0.3 in every chain. The report also gives slope MCSE. These criteria follow [Stan's diagnostic guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html). This project's additional conservative gate requires no tree-depth saturation. All assess computation, not substantive adequacy.

If a non-smoke fit fails its diagnostic gate, the pipeline retains the original diagnostics and makes one retry with 10,000 iterations, 5,000 warmup iterations, and `adapt_delta = 0.99`. The maximum tree depth increases to 14 only if the original fit had tree-depth hits. If the retry still fails, a targeted third attempt raises maximum tree depth to 14 when saturated, and extends to 20,000 iterations with 5,000 warmup when R-hat or ESS still fails; adapt_delta remains 0.99. Any remaining failure stops the run as unverified. These computational extensions do not change priors, data, or model structure and do not use predictive performance. No model is silently removed. The per-fit metadata records all attempts and the settings actually used for the accepted output.

Posterior predictive checks condition on observed patients and, for H/C, their posterior random effects. They examine the FVC median, 5th/95th percentiles and standard deviation; the across-patient standard deviations of patient means and ordinary-least-squares slopes; the mean patient slope; pooled adjacent-visit residual correlation; and pooled time slopes of signed and absolute residuals. Patient slopes require at least two distinct times. Residual pairs join successive visits within each patient after sorting by time, and do not assume equal intervals. Patient and residual slope statistics in these checks are in mL/year, as indicated in their names.

Observed distributional statistics are fixed, while observed residual discrepancies depend on each draw's fitted mean. For residual checks, the public `observed` entry is the median of those draw-specific discrepancies. The predictive tail fraction compares replicated and observed discrepancies within the same draw before averaging. Aggregate predictive checks do not establish that every individual's trajectory is well represented.

The report presents the primary inference alongside both sensitivity analyses, without selecting whichever result appears most favorable. Covariate associations, especially those involving small subgroups, are exploratory. Unequal follow-up, informative dropout, cohort selection, measurement variability, and linear/Gaussian assumptions can affect the results. Age adjustment and smoking indicators do not turn these associations into causal effects.
