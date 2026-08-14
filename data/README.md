# Input data

The CSV files the analysis reads are not distributed in this repository. They contain the
aggregate output of queries run against the InEK DatenBrowser, which is a public query
interface rather than a downloadable dataset; the outputs are reproducible by anyone from the
query specification in `../protocol/query_grid.csv`.

Expected files, with the column specification in `../analysis/schema/input_schema.md`:

- `arm_totals.csv`
- `arm_distributions.csv`
- `hospital_tables.csv`
- `procedure_family_counts.csv`
- `secondary_diagnosis_totals.csv`
- `secondary_diagnosis_distributions.csv`
- `cohort_overlap_totals.csv`
- `diagnosis_specificity.csv`
- `secondary_diagnosis_panel.csv`
- `comorbidity_exact.csv`
- `discharge_destination.csv`

The first three are required; the analysis skips the sections that depend on any of the others
if the file is absent.

Once they are in place, `python3 ../analysis/01_analysis.py` regenerates everything in
`results/`, and `04_verify.py` checks the transcription and the derived quantities.
