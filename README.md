# Bayesian Modeling of Pulmonary Fibrosis Trajectories

**Junguang Jia**

A reproducible, exploratory study of longitudinal forced vital capacity (FVC), implemented in R and Stan. The analysis estimates average change and between-patient variation, examines associations with age, sex, and smoking status, and evaluates predictions for patients with no observed FVC history.

Three prespecified Gaussian models connect the scientific and predictive questions:

| Model | Structure | Purpose |
|:--|:--|:--|
| P | Pooled intercept and time slope | Descriptive benchmark |
| H | Correlated patient-specific intercepts and slopes | Primary inference on longitudinal change |
| C | H plus age, sex, smoking status, and their time interactions | Exploratory adjusted associations |

Five-fold validation holds out entire patients. Predictive scores integrate out the held-out patient's random effects. This is internal validation within the supplied cohort; it is not external validation or evaluation of a clinical decision tool.

## Research release

- [Research report (PDF)](results/research/report.pdf) and [standalone LaTeX source](results/research/report.tex).
- [Computational validation report (PDF)](results/research/validation.pdf) and [LaTeX source](results/research/validation.tex).
- [Run manifest](results/research/manifest.json), [sampling diagnostics](results/research/diagnostics.csv), and [patient-level validation aggregates](results/research/cv_metrics.csv).

The retained research run contains 176 patients and 1,549 measurements, consolidated to 1,542 patient-weeks. All 24 model configurations passed the recorded computational gates after 34 sampling attempts. The primary H mean slope is **-4.25 mL/week** (95% credible interval **[-5.12, -3.40]**). This estimate is conditional on the model: posterior checks reveal remaining residual serial dependence. Passing sampling diagnostics does not resolve that limitation or establish external validity.

The report can be regenerated from the committed aggregate summaries with `Rscript scripts/render.R --run research`. Re-fitting requires an independently authorized local data copy.

## Start with synthetic data

Requirements: R, a working C++ toolchain for RStan, and a LaTeX distribution with `pdflatex` for PDF reports. Reports use standard LaTeX packages and contain their figures directly in the source. Use `--tex-only` to generate a standalone source for another LaTeX editor. Install the pinned project dependencies, run unit checks, and generate a short synthetic demonstration:

```sh
Rscript scripts/setup.R
Rscript tests/run_tests.R
Rscript scripts/run.R --mode demo --smoke
Rscript scripts/render.R --run demo
```

The smoke run checks compilation and execution with reduced sampling. It does not establish convergence or support scientific conclusions. Omit `--smoke` for the full synthetic demonstration.

## Reproduce the research analysis

Obtain an authorized local CSV and review the [data requirements and provenance note](docs/data.md). Then run:

```sh
Rscript scripts/run.R --mode research --data /path/to/authorized-data.csv
Rscript scripts/render.R --run <run-id>
```

The pipeline reports its run identifier. Public summaries are written to `results/<run-id>/`, while private data, fits, logs, and patient-level calculations stay outside the public artifact set. Rendering uses the selected run's summaries, never another project's cache. The run manifest records data and configuration hashes, package versions, actual sampling settings, and analysis status.

Research defaults are four chains, 6,000 iterations per chain including 3,000 warmup iterations, seed 4224, `adapt_delta = 0.95`, and `max_treedepth = 12`. The prespecified research analysis comprises 24 model configurations: three full-data models, 15 cross-validation fits, four prior-sensitivity fits, and two duplicate-measurement sensitivity fits. Failed sampling gates trigger recorded computational extensions; the number of sampling attempts can therefore exceed 24. Runtime depends on the local toolchain and hardware.

A separate seeded parameter-recovery check uses 100 synthetic participants:

```sh
Rscript scripts/run.R --mode recovery
Rscript scripts/render.R --run recovery
```

It checks the adjusted overall slope against a known generating value with a prespecified absolute-error tolerance of 1.5 mL/week. This single experiment is not a calibration study.

## Methods and outputs

- [Methods](docs/methods.md): estimands, priors, prediction targets, diagnostics, and sensitivity analyses.
- [Data documentation](docs/data.md): required fields, preprocessing, access, and publication boundaries.
- [Report generator](scripts/render.R): run-specific standalone LaTeX source and PDF, generated from verified aggregate summaries.
- `results/`: public tables and rendered reports for retained runs, when available.

Each report distinguishes synthetic demonstrations from research runs and reports computational failures explicitly. Inference requires the prespecified sampling diagnostics to pass; passing those diagnostics alone does not establish model adequacy. Predictive checks and sensitivity analyses remain part of interpretation.

## Data access and privacy

No real participant-level CSV, participant identifier, full fitted object, or individual prediction is distributed. The supplied data schema resembles the tabular component of the [OSIC Pulmonary Fibrosis Progression competition](https://www.kaggle.com/competitions/osic-pulmonary-fibrosis-progression/data), but exact provenance of the local input has not been independently established. Obtain data directly from the authorized provider and follow the applicable [competition rules](https://www.kaggle.com/competitions/osic-pulmonary-fibrosis-progression/rules). The code license does not grant rights to third-party data.

## Interpretation and reuse

The project is observational and exploratory. Covariate coefficients describe model-based associations, not causal effects. Small subgroups, linear trajectories, informative follow-up, and Gaussian residual assumptions limit the conclusions. Population-mean credible intervals and new-patient predictive intervals answer different questions and are labelled separately.

Code and original project documentation are released under the [MIT License](LICENSE). Use the metadata in [CITATION.cff](CITATION.cff) to cite the software. Cite the original data provider separately when using their data.
