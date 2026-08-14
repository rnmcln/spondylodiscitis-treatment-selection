# Operative treatment and in-hospital mortality in pyogenic spondylodiscitis

Query protocol, analysis code and derived results for a national analysis of 61,134 German
discharge episodes with a principal diagnosis of pyogenic spondylodiscitis (ICD-10-GM
M46.2–M46.5), data years 2019 to 2024.

## Study summary

**Question.** Who undergoes surgery for pyogenic spondylodiscitis across an entire national
hospital system, how do the operated and non-operated populations differ, and what does that
selection mean for the observed difference in in-hospital mortality between them?

**Data.** The InEK DatenBrowser (<https://datenbrowser.inek.org>) is a public query interface
to the §21 KHEntgG dataset, the mandatory annual submission of case-level billing records by
German hospitals reimbursed under the aG-DRG system. It returns aggregate counts and
distributions, not case-level records. The unit of observation is the acute inpatient
discharge episode; episodes cannot be linked to patients, and cells with fewer than four
episodes are suppressed.

**Design.** Because the interface returns aggregates, every quantity is the difference of two
queries. Four queries per data year define the 2 × 2 table (cohort, cohort with an operative
code, cohort with in-hospital death, cohort with both); non-operative quantities follow by
subtraction. Ten further queries per year cover complexity and resource-use markers, two
alternative exposure definitions and an external anchor, giving 84 queries in total. Two
further series were added: one query per operative OPS family per year (54), the four-query
structure repeated on the secondary-diagnosis population (24) and on the overlap between the
two (24), the distributions of the secondary-diagnosis cohort (24), a split of the cohort by
diagnosis specificity (24), a read of the secondary-diagnosis register in each of the four
query states (24), filter queries for five diagnosis rubrics by arm (72) and a discharge
destination series (48) and a debridement-directed exposure definition (12).
`protocol/` contains the full specification.

**Exposure.** An episode is operative if at least one code from nine prespecified OPS families
covering spinal decompression, debridement, reconstruction and instrumentation was recorded
during the admission, each family expanded to all terminal subcodes valid in that OPS year
(174 to 218 codes per year). The exposure is therefore "operated at any point during the
episode": it carries no treatment intention and no procedure timing.

**Main findings.**

| Quantity | Value |
|---|---|
| Discharge episodes, 2019–2024 | 61,134 |
| Episodes with an operative code | 26,130 (42.7%) |
| Operative share, age 50–54 → ≥80 years | 52.6% → 33.2% |
| Operative share, <300 beds → ≥1000 beds | 27.9% → 65.2% |
| PCCL 4–6, operative vs non-operative | 49.3% vs 13.7% |
| Intensive care, operative vs non-operative | 50.5% vs 10.9% |
| In-hospital mortality, operative vs non-operative | 7.05% vs 6.30% |
| Crude risk ratio | 1.12 (95% CI 1.05–1.19) |
| Age-adjusted risk ratio (Mantel-Haenszel) | 1.31 (95% CI 1.24–1.39) |
| Intraspinal abscess (G06.1), operative vs non-operative | 22.6% vs 5.8% |
| Sepsis (A41), operative vs non-operative | 9.9% vs 5.5% |
| Risk ratio stratified by year and abscess | 1.05 (0.99–1.12) |
| Risk ratio stratified by year and sepsis | 0.93 (0.88–0.98) |
| Transfer to another hospital, operative vs non-operative | 18.5% vs 25.1% |
| Directly age-standardised mortality difference | 1.84 percentage points (1.43–2.25) |
| Secondary-diagnosis-only cohort (39,239 episodes) | 0.94 (0.89–1.00) |
| Secondary-diagnosis cohort including the overlap (60,244) | 0.84 (0.80–0.88); age-adjusted 0.93 (0.88–0.97) |
| Union of both cohorts (100,373 episodes) | 0.93 (0.89–0.96) |
| Cohort carrying the specific codes M46.2/M46.3 | 4524 (7.4%), mortality 6.87% vs 6.60% |
| E-value for the age-adjusted estimate | 1.95 |

**Interpretation.** Treatment allocation is structured along several axes at once, and they
do not point the same way.

- *Age.* Operated patients are younger, and younger age predicts lower mortality, so adjusting
  for age moves the association further from the null rather than towards it.
- *Acute coded conditions.* Intraspinal abscess, sepsis, septic shock and paraparesis are all
  enriched in the operative arm. Stratifying by any of them moves the association in the
  opposite direction to age adjustment.
- *Chronic comorbidity.* Heart failure and chronic kidney disease run the other way, being
  commoner in the non-operative arm; diabetes is equally frequent in both.
- *Cohort definition.* In episodes carrying M46.2–M46.5 only as a secondary diagnosis, where
  the admission is driven by another condition, mortality is twice as high overall and the
  association is reversed.

Three properties of the data bound what any of this can support. The interface returns
marginal tables, so no two stratifications can be combined and no estimate here is mutually
adjusted. Secondary diagnoses describe the whole admission, so a coded condition may follow an
operation as readily as precede it, and none of them is a baseline risk-adjustment covariate.
In-hospital death is censored by transfer to another hospital, which ends a quarter of
non-operative and a fifth of operative episodes, and post-transfer outcomes are not observable.

The estimates therefore describe the association between operative coding and in-hospital
death, and how strongly that association depends on which dimension of selection and which
episode population is represented. They do not estimate an effect of surgery and are not
evidence of benefit or harm.

## Repository layout

```
protocol/
  exposure_codebook.csv     OPS families, role in each exposure definition, rationale
  query_grid.csv            every query run, with its entered parameters
  extraction_protocol.md    source, cohort, exposure, verification rules, suppression
analysis/
  01_analysis.py            association measures, standardisation, sensitivity analyses
  02_figures.py             figures 1-3 and supplementary figures S1-S2
  03_tables.py              manuscript tables
  04_verify.py              consistency checks on inputs and results
  schema/input_schema.md    column specification of the input files
results/                    derived summary results and tables, including per-OPS-family
                            counts and the cohort-definition sensitivity
figures/                    figures 1-3 and S1-S2 as vector PDF and 300 dpi PNG
data/                       input files (not distributed; see below)
```

## Data availability

The DatenBrowser query outputs are not redistributed here. `data/` therefore ships empty
apart from its README. To reproduce the analysis, run the queries in
`protocol/query_grid.csv` against the public interface and transcribe the results into the
input files described in `analysis/schema/input_schema.md`. The derived summary results
in `results/` are the quantities reported in the paper and are provided in full.

## Reproducing

```bash
python3 analysis/01_analysis.py    # writes results/
python3 analysis/02_figures.py     # writes figures/
python3 analysis/03_tables.py      # writes the tables
python3 analysis/04_verify.py      # 175 consistency checks; non-zero exit on failure
```

Python 3.11. `01_analysis.py`, `03_tables.py` and `04_verify.py` use the standard library
only; `02_figures.py` requires matplotlib (see `requirements.txt`).

## Methods notes

- **Age strata.** The interface publishes age distributions as percentages to two decimal
  places. Stratum counts are reconstructed from those percentages and rescaled so that each
  stratum set sums exactly to its query total. Two sensitivity analyses bound the effect:
  resampling every percentage within its rounding interval (2000 draws) moves the pooled
  age-adjusted risk ratio between 1.3129 and 1.3132.
- **Uncertainty.** Confidence intervals treat the observed episodes as realisations of an
  underlying care process. They do not account for repeated admissions by the same patient,
  which the data cannot identify.
- **PCCL** is derived by the annual aG-DRG grouper and is not comparable across data years; it
  is used descriptively and not for adjustment.
- **Suppression.** Suppressed residual cells (566 episodes in Q1, 73 in Q2, 10 in Q3, 32 in
  Q4) are excluded from the hospital-structure analysis and reported rather than imputed. The
  secondary-diagnosis register lists codes down to five episodes; a query returning too few
  cases is not displayed at all, which is why the hospice discharge counts are reported for
  the whole cohort only.
- **Coded conditions are not baseline covariates.** Secondary diagnoses describe the whole
  admission, so a condition may precede an operation, prompt it or follow it. Stratified
  estimates using them are alternatives to age adjustment, not refinements of it, and they
  cannot be combined with it.
- **Cohort overlap.** 21,005 episodes carry one M46.2–M46.5 code as principal diagnosis and
  another as secondary, so they belong to both cohorts. The overlap was measured directly and
  removed, which makes the secondary-only cohort (39,239 episodes) and the union (100,373)
  exact.
- **Rubric counts.** A count for an ICD rubric with several subcodes cannot be obtained by
  summing register entries, because an episode may carry more than one subcode. Measured
  against filter queries the register sums were 1.5% (N18), 2.3% (E11), 4.6% (A41), 32.2%
  (I50) and 68.1% (G82) too high. Rubric counts therefore come from filter queries, and only
  codes without subcodes are read from the register.

## Licence

Code: MIT (`LICENSE`). Derived results and documentation: CC BY 4.0.
