# Data, access, and preprocessing

## Access and provenance

Real participant-level data are not included in this repository. The field names and longitudinal structure of the available input resemble the tabular data in the [OSIC Pulmonary Fibrosis Progression competition](https://www.kaggle.com/competitions/osic-pulmonary-fibrosis-progression/data). Exact identity with an original provider file has not been independently established; the project must not be cited as an independently verified redistribution of that dataset.

For reproduction, obtain an authorized copy directly from the provider, review the applicable [competition rules](https://www.kaggle.com/competitions/osic-pulmonary-fibrosis-progression/rules), and supply its local path to the analysis command. Access conditions and permissions belong to the provider. The MIT license covers this project's code and original documentation, not third-party data.

The run manifest includes a SHA-256 digest of the actual input, permitting local comparisons without distributing it. The study report derives all cohort counts from the selected run. It does not assume that another provider export has the same number of records or participants.

## Input dictionary

The CSV uses the following case-sensitive fields:

| Field | Expected contents | Analysis use |
|:--|:--|:--|
| `Patient` | Nonmissing patient identifier | Private grouping key; never published |
| `Weeks` | Finite numeric time in weeks | `time = Weeks / 52` |
| `FVC` | Positive finite forced vital capacity in mL | `outcome = FVC / 1000`, in litres |
| `Age` | Positive finite numeric age in years | `(Age - 65) / 10` |
| `Sex` | `Female` or `Male` | Patient-level covariate |
| `SmokingStatus` | `Never smoked`, `Ex-smoker`, or `Currently smokes` | Patient-level covariate |
| `Percent` | Optional numeric percent-predicted FVC, when supplied | Excluded from every model |

Age, sex, and smoking status must be constant within each patient. Missing or invalid model inputs are errors, rather than a silent complete-case analysis. Exact duplicate analysis records are rejected for source review; repeated patient-week measurements with different FVC values are handled below. Extra columns are not predictors. The input validation establishes that the requested analysis can run; it does not authenticate the origin of the data.

The original time origin is retained. Negative weeks are permitted. Time zero must not be interpreted as diagnosis, symptom onset, or study entry without additional source documentation.

## Duplicate measurements

The main analysis averages FVC values for rows sharing the same patient and week. Each resulting patient-week contributes one outcome with equal residual variance in the model. Counts of input rows, unique patients, unique patient-weeks, and removed duplicate rows are calculated before fitting and retained as aggregate summaries.

A sensitivity analysis refits H and C using every original row. This evaluates the consequence of the duplicate-handling choice; it does not prove that repeated same-week measurements are conditionally independent.

## Publication boundaries

Public artifacts contain aggregate cohort counts, population-level posterior summaries, computational diagnostics, aggregate predictive metrics and checks, population-mean curves, and a run manifest. They exclude raw rows, patient identifiers or their hashes, patient-specific random effects, per-patient predictions or scores, and full posterior fit objects. The report contains no original patient trajectories or row-level tables.

The synthetic demonstration is generated independently of real participant records. It is suitable for testing the workflow and does not reproduce an individual's measurements. Even when its software checks pass, a synthetic run is not evidence about pulmonary fibrosis in the research cohort.

## Included synthetic example

[`data/synthetic/example.csv`](../data/synthetic/example.csv) contains 15 generated participants and 120 rows from `simulate_data(15L, 4224L)`. Every grouping key begins with `SYNTHETIC-`. This file is generated from specified distributions without reading any real data and is released under the project's MIT license. The demo command regenerates its own synthetic input, so it never substitutes this example for research observations.
