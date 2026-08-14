# Input file schema

The analysis reads the CSV files below from `data/`. They hold the aggregate output of the
queries in `protocol/query_grid.csv`. The DatenBrowser output is not redistributed in this
repository; rebuild these files by running the query grid against the public interface.

## `data/arm_totals.csv`

One row per data year.

| Column | Content |
|---|---|
| `year` | data year, 2019–2024 |
| `N_total` | Q1: episodes with the principal diagnosis |
| `N_operative` | Q2: episodes with at least one core operative code |
| `N_nonoperative` | Q1 − Q2 |
| `deaths_total` | Q3: episodes ending in in-hospital death |
| `deaths_operative` | Q4 |
| `deaths_nonoperative` | Q3 − Q4 |
| `N_operative_wide`, `deaths_operative_wide` | S1, S2: wide exposure definition |
| `N_instrumented`, `deaths_instrumented` | S3, S4: instrumented exposure definition |
| `ICU_operative`, `ICU_total` | M1, M2: episodes with intensive care |
| `vent96_operative`, `vent96_total` | M3, M4: episodes with ventilation from 96 hours |
| `reop_operative` | M5: operative episodes carrying OPS 5-983 |

## `data/arm_distributions.csv`

One row per data year and query state. `state` takes the values `cohort` (Q1), `operative`
(Q2), `deaths_total` (Q3) and `deaths_operative` (Q4).

| Column | Content |
|---|---|
| `year`, `state` | key |
| `fallzahl` | episode count of that query, used to check the transcription |
| `male_pct` | percentage male |
| `vwd_mean`, `vwd_sd` | mean and SD of length of stay, days |
| `kurz_pct`, `normal_pct`, `lang_pct` | short-stay, normal-stay and long-stay percentages |
| `pccl0` … `pccl6` | PCCL distribution, percentages |
| `age_u30`, `age3039`, `age4049`, `age5054`, `age5559`, `age6064`, `age6574`, `age7579`, `age80p` | age distribution, percentages |

The interface prints the age distribution in two blocks and sixteen bands. The eight bands
below 30 years (`< 28 Tage`, `28 Tage - 1. Jahr`, `1-2 Jahre`, `3-5 Jahre`, `6-9 Jahre`,
`10-15 Jahre`, `16-17 Jahre`, `18-29 Jahre`) are summed into `age_u30`; the remaining eight
are transcribed directly. Each row must therefore sum to 100 within the rounding of sixteen
two-decimal values (about 0.08 percentage points); `04_verify.py` enforces this.

All percentages are transcribed exactly as published, to two decimal places. Each
distribution must sum to 100 within transcription rounding.

## `data/hospital_tables.csv`

The bed-size class × ownership cross-tabulation, one row per cell per query per data year.

| Column | Content |
|---|---|
| `year` | data year |
| `query` | `Q1`, `Q2`, `Q3` or `Q4` |
| `bettengroessenklasse` | bed-size class label as printed by the interface |
| `traegerschaft` | ownership label as printed (`Oeffentlich`, `Freigemeinnuetzig`, `Privat`) |
| `faelle` | episode count in the cell |
| `suppressed_rest` | `True` for the residual row that pools suppressed cells |

For each year and query the `faelle` column must sum exactly to the corresponding total in
`arm_totals.csv`.

## `data/procedure_family_counts.csv`

One row per data year and operative OPS family (54 rows).

| Column | Content |
|---|---|
| `year` | data year |
| `family` | OPS family, e.g. `5-836` |
| `codes` | terminal codes selected for that family in that year |
| `episodes` | episodes carrying at least one code from that family |

Families overlap; the counts sum to more than the operative arm.

## `data/secondary_diagnosis_distributions.csv`

The same four states as `arm_distributions.csv`, for the secondary-diagnosis cohort: one row
per data year and state (`cohort`, `operative`, `deaths_total`, `deaths_operative`), same
columns. Used to age-adjust the cohort-definition comparison.

## `data/secondary_diagnosis_totals.csv`

One row per data year, for the cohort-definition sensitivity analysis.

| Column | Content |
|---|---|
| `year` | data year |
| `N_total` | episodes with M46.2–M46.5 as a secondary diagnosis |
| `N_operative` | of those, episodes with at least one core operative code |
| `deaths_total`, `deaths_operative` | corresponding in-hospital deaths |

## `data/cohort_overlap_totals.csv`

Episodes carrying an M46.2–M46.5 code as principal **and** another as secondary diagnosis,
one row per data year, same four count columns as `secondary_diagnosis_totals.csv`. Used to
compute the exact union of the two cohorts.

## `data/diagnosis_specificity.csv`

One row per data year, splitting the primary cohort by diagnosis specificity.

| Column | Content |
|---|---|
| `year` | data year |
| `N_specific`, `deaths_specific` | episodes and deaths with M46.2 or M46.3 as principal diagnosis |
| `N_unspecific`, `deaths_unspecific` | episodes and deaths with M46.4 or M46.5 |

Each year's two parts must sum exactly to the cohort and death totals in `arm_totals.csv`;
`04_verify.py` enforces this, and the extraction used the same rule as its acceptance test.

## `data/secondary_diagnosis_panel.csv`

The secondary-diagnosis register of each of the four queries, one row per data year and query
state (24 rows). The register lists every ICD-10-GM code carried by at least five episodes of
that query, with the number of episodes carrying it; the columns below hold the subset used in
the analysis.

| Column | Content |
|---|---|
| `year`, `q` | data year and query state (`Q1`–`Q4`) |
| `n` | case count implied by the register (count divided by its published percentage) |
| `nenn` | total code entries (`Nennungen`) across all listed codes and the suppressed residual |
| `meanDx` | `nenn` divided by `n`, the mean number of coded secondary diagnosis entries per episode |
| `restF` | episodes carrying at least one code that fell below the suppression floor |
| `faelle` | total episode counts (`Fälle`) summed over all listed codes, plus `restFaelle` |
| `restFaelle` | the suppressed residual row's episode count, equal to `restF` |
| `codes` | number of codes listed in that register read |
| remaining columns | episodes carrying a code in the named rubric, summed over its subcodes |

`nenn` counts entries, not distinct codes: an episode contributes more than one entry for a
code recorded in more than one department segment (in the 2024 cohort, A41.51 appears as 71
episodes and 114 entries). `faelle` divided by `n` is a distinct-code count, but a lower
bound, because the suppressed residual row reports episodes carrying at least one suppressed
code rather than code-episode pairs.

The summed rubric columns are **not** episode counts for the rubric: an episode may carry more
than one subcode of the same rubric, and codes below the suppression floor are absent. Only
columns for codes with no subcodes (`G06.1`, `R57.2`, `G06.0`, `G06.2`, `D62`, `N19`, `Z92.1`)
are exact episode counts. Rubric-level counts come from `comorbidity_exact.csv`.

## `data/comorbidity_exact.csv`

Filter-query counts for diagnosis rubrics, one row per data year, arm, rubric and outcome.

| Column | Content |
|---|---|
| `year` | data year |
| `arm` | `cohort` (principal diagnosis only) or `operative` (principal diagnosis + operative codes) |
| `group` | ICD rubric entered in the secondary-diagnosis field, e.g. `A41` |
| `outcome` | `all`, or `death` for the additional discharge-reason filter |
| `n_codes` | terminal codes selected for that rubric in that data year |
| `n` | episodes returned |

Non-operative counts follow as `cohort − operative`. Sepsis (`A41`) carries the `death` rows,
which permit a 2 × 2 × 2 table by arm, sepsis and death within each data year.

## `data/discharge_destination.csv`

One row per data year, arm and destination group.

| Column | Content |
|---|---|
| `year` | data year |
| `arm` | `cohort` or `operative` |
| `destination` | discharge-reason codes entered, `06+08` (transfer to another acute hospital), `09` (rehabilitation), `10` (nursing facility), `11` (hospice) |
| `n` | episodes returned, or `SUP` when the interface reported too few cases to display |

## `data/exposure_debridement.csv`

One row per data year, for the debridement-directed exposure sensitivity.

| Column | Content |
|---|---|
| `year` | data year |
| `codes` | terminal codes selected across the seven families in that year |
| `N_debridement` | episodes with at least one code from 5-033, 5-830 to 5-832 or 5-835 to 5-837 |
| `deaths_debridement` | of those, episodes ending in in-hospital death |

The definition excludes the generic family 5-839 and the dominant 5-83b. Its comparator is the
non-operative arm as defined by the core exposure, so episodes carrying only an excluded
operative code fall outside both arms of that comparison.
