# Extraction protocol

## Source

InEK DatenBrowser, <https://datenbrowser.inek.org>. The DatenBrowser is a public query
interface to the §21 KHEntgG dataset, the mandatory annual submission of case-level billing
records by German hospitals reimbursed under the aG-DRG system. It returns aggregate counts
and distributions for a filter combination; it does not return case-level records and it
does not offer a bulk download.

Data years 2019 to 2024. Unit of observation: the acute inpatient discharge episode.

## Cohort

Principal diagnosis (Hauptdiagnose) ICD-10-GM M46.2, M46.3, M46.4 and M46.5, entered as all
40 terminal five-digit codes. Secondary diagnoses are not used.

## Exposure

An episode is operative if at least one code from the core OPS families in
`exposure_codebook.csv` was recorded during the admission. Each family is expanded to every
terminal subcode valid in the OPS catalogue of the respective data year, giving 174 codes
(2019), 180 (2020), 218 (2021), 218 (2022), 212 (2023) and 214 (2024). Selection within the
procedure field is disjunctive; the procedure field combines conjunctively with the
diagnosis and discharge-reason fields.

Two additional definitions are used for sensitivity analysis:

- **wide**: core plus the lumbar access family 5-032 (199 codes in 2019).
- **instrumented**: 5-835, 5-836, 5-837 and 5-83b only (87 codes in 2019), which excludes
  the generic families 5-839 and 5-033.

## Outcome

In-hospital death, selected on discharge reason (Entlassungsgrund) "07 Tod".

## Query grid

Fourteen queries per data year, 84 in the primary grid, with the later series listed
alongside them in `query_grid.csv`.

| Query | Filter | Purpose |
|---|---|---|
| Q0 | cohort, read from the Info register | external anchor against Destatis |
| Q1 | cohort | cohort size and distributions |
| Q2 | cohort + core exposure | operative arm |
| Q3 | cohort + death | total deaths |
| Q4 | cohort + core exposure + death | operative deaths |
| M1, M2 | cohort ± exposure, intensive care | resource-use marker |
| M3, M4 | cohort ± exposure, ventilation from 96 h | resource-use marker |
| M5 | cohort + exposure + procedure 5-983 | reintervention marker |
| S1, S2 | wide exposure, with and without death | exposure sensitivity |
| S3, S4 | instrumented exposure, with and without death | exposure sensitivity |

Non-operative quantities are obtained by subtraction: Q1 − Q2 and Q3 − Q4.

Two further series were run after the primary grid:

- **F1–F9**, one query per operative OPS family per data year (54 queries), giving the number
  of episodes identified by each family alone. Selection within the procedure field is
  disjunctive, so these counts are not recoverable from the pooled query; families overlap and
  their counts sum to more than the operative arm.
- **X1–X4**, the four-query structure repeated with M46.2–M46.5 entered as *secondary*
  diagnosis and no principal-diagnosis restriction (24 queries), for the cohort-definition
  sensitivity analysis. Totals and full distributions were captured, so that comparison is
  age-adjusted by the same method as the primary cohort.
- **O1–O4**, the same four-query structure with M46.2–M46.5 entered in *both* the principal
  and the secondary diagnosis field (24 queries). Because every episode counted by this series
  belongs to the secondary-diagnosis cohort, subtracting it yields the secondary-only cohort
  (episodes in which spondylodiscitis is never the principal diagnosis) exactly, and it also
  makes the union of the two cohorts exact.
- **D1–D2**, the cohort and death counts split into the specific codes (M46.2, M46.3) and the
  unspecified codes (M46.4, M46.5), 24 queries. Each split was accepted only when its two
  parts summed exactly to the corresponding annual total, which also served as the
  verification rule for the selection.

- **N1–N4**, the four-query structure re-run and read from the **Nebendiagnosen** register
  rather than the Info register (24 queries). The register lists every ICD-10-GM secondary
  diagnosis carried by at least five episodes of that query, with the number of episodes
  carrying it, so arm-specific and death-specific counts follow by the same subtraction. Each
  read was accepted only when the case count implied by the first row (count divided by its
  published percentage) matched the query total. The entry percentages were summed as a
  secondary check; across the six cohort reads they lay between 99.06 and 99.82, which is
  what rounding across about 1900 values printed to two decimal places produces.
- **C-series**, filter queries entering an ICD rubric in the secondary-diagnosis field, for
  A41, G82, E11, I50 and N18, in the cohort and the operative arm of every data year, with the
  additional death filter for A41 (72 queries). These are needed because a rubric count cannot
  be obtained by summing register entries: an episode may carry more than one subcode of the
  same rubric. Measured against these queries the register sums were 1.5% (N18), 2.3% (E11),
  4.6% (A41), 32.2% (I50) and 68.1% (G82) too high.
- **B1–B2**, a debridement-directed exposure definition (5-033, 5-830 to 5-832, 5-835 to
  5-837), with and without the death filter, for every data year (12 queries). It drops the
  generic family 5-839 and the dominant 5-83b, which appears in 85.7% of operative episodes,
  so that the result can be checked against a definition built only from procedures aimed at
  the infected segment.
- **Discharge series**, the cohort and the operative arm of every data year filtered on
  discharge reason 06+08 (transfer to another acute hospital), 09 (rehabilitation), 10
  (nursing facility) and 11 (hospice), 48 queries. The operative hospice count fell below the
  display threshold of the interface in every data year and is reported for the whole cohort
  only.

Each family query was verified before it was run: every code held in the procedure field had
to carry the family prefix, and the per-year terminal-code counts had to sum to the core code
count for that year (174, 180, 218, 218, 212 and 214). Three family selections failed this
check on the first attempt and were rebuilt.

## Verification rules applied during extraction

1. Every cross-tabulated distribution returned by the interface (age, sex, PCCL, length of
   stay, bed-size class × ownership) must sum exactly to the case count of its own query.
   All 24 bed-size × ownership tables reconciled exactly.
2. Each Q3 result is reproduced after the procedure filter has been removed and rebuilt from
   scratch, to detect residual filter state in the interface. This rule identified one
   contaminated marker value (2023, intensive care), which was re-extracted.
3. Annual cohort totals are compared with the hospital diagnosis statistics of the Federal
   Statistical Office (Destatis) for the same diagnosis group and calendar year.
4. Every transcribed distribution must sum to 100 within the rounding of the published
   values. The age distribution is printed in two blocks of eight bands; all sixteen must be
   captured, with the eight bands below 30 years summed into a single stratum.
5. The whole filter state was read back and compared with the intended query immediately
   before each run, including the fields that were supposed to be empty. Residual state from a
   previous query is the failure mode that produces plausible but wrong counts, and it is not
   otherwise visible.
6. Every code selection was verified by prefix membership rather than by count: each code held
   in the field had to carry the intended prefix before the query was run.

## Suppression

Cells with fewer than four episodes are suppressed by the interface and pooled into a
residual ("Rest") row. The secondary-diagnosis register lists codes down to five episodes, and
a query whose result is smaller than the display threshold returns no table at all rather than
a suppressed cell. Across the study period the suppressed residuals amounted to 566
episodes in Q1, 73 in Q2, 10 in Q3 and 32 in Q4. They are excluded from the hospital-structure
analysis and reported rather than imputed.

## Reproducing the input files

The DatenBrowser output is not redistributed here. To rebuild the input files
described in `../analysis/schema/input_schema.md`, run each query in `query_grid.csv`
against the live interface and transcribe the returned totals and distributions into the
documented columns.
