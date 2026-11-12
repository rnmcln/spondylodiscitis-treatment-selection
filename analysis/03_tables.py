#!/usr/bin/env python3
"""Build the manuscript tables from the analysis outputs.

Writes results/table1.csv, results/table2.csv, results/table3.csv and
results/tables.json. All rounding is half-up so that the printed values are
reproducible from the four-decimal analysis outputs.
"""
import csv, json, os
from decimal import Decimal, ROUND_HALF_UP

RES = os.environ.get('SPONDYLO_RESULTS', os.path.join(os.path.dirname(__file__), '..', 'results'))
rd = lambda n: list(csv.DictReader(open(os.path.join(RES, n))))


def r(x, dp=2):
    q = Decimal(1).scaleb(-dp)
    return str(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP))


def thousands(n):
    n = int(round(float(n)))
    return f'{n:,}' if n >= 10000 else str(n)


arm = rd('arm_comparison.csv')
rr = {x['year']: x for x in rd('risk_ratios.csv')}
ann = {x['year']: x for x in rd('annual_summary.csv')}
age = rd('age_specific.csv')
hosp = rd('hospital_structure.csv')
std = {x['quantity']: x for x in rd('standardised_mortality.csv')}
sens = rd('exposure_sensitivity.csv')
smr = {x['year']: x for x in rd('smr.csv')}

# ---------------------------------------------------------------- Table 1
t1 = [['Characteristic', 'Operative (n = 26,130)', 'Nonoperative (n = 35,004)']]
for row in arm:
    lab, o, c = row['characteristic'], row['operative'], row['non_operative']
    if lab == 'Admissions, n':
        continue
    if lab == 'In-hospital deaths, n':
        o, c = thousands(o), thousands(c)
    elif lab == 'In-hospital mortality, %':
        o, c = r(o), r(c)
    lab = lab.replace('Age >=80 years', 'Age ≥80 years').replace('PCCL 4-6', 'PCCL 4–6')
    if lab.startswith('Age 65-79'):
        lab = lab.replace('65-79', '65–79')
    t1.append([lab, str(o), str(c)])

# ---------------------------------------------------------------- Table 2
t2 = [['Data year', 'Operative, deaths/n (%)', 'Nonoperative, deaths/n (%)',
       'Crude RR (95% CI)', 'Age-adjusted RR (95% CI)']]
for y in [x['year'] for x in rd('annual_summary.csv')]:
    a = ann[y]
    N, O = int(a['episodes']), int(a['operative'])
    D, DO = int(a['deaths']), int(a['deaths_operative'])
    q = rr[y]
    t2.append([y,
               f"{DO}/{thousands(O)} ({r(a['mort_operative_pct'])})",
               f"{D-DO}/{thousands(N-O)} ({r(a['mort_nonoperative_pct'])})",
               f"{r(q['crude_rr'])} ({r(q['crude_lo'])}–{r(q['crude_hi'])})",
               f"{r(q['mh_rr'])} ({r(q['mh_lo'])}–{r(q['mh_hi'])})"])
tot_O = sum(int(ann[y]['operative']) for y in ann)
tot_N = sum(int(ann[y]['episodes']) for y in ann)
tot_D = sum(int(ann[y]['deaths']) for y in ann)
tot_DO = sum(int(ann[y]['deaths_operative']) for y in ann)
q = rr['pooled']
t2.append(['Pooled',
           f"{thousands(tot_DO)}/{thousands(tot_O)} ({r(100*tot_DO/tot_O)})",
           f"{thousands(tot_D-tot_DO)}/{thousands(tot_N-tot_O)} ({r(100*(tot_D-tot_DO)/(tot_N-tot_O))})",
           f"{r(q['crude_rr'])} ({r(q['crude_lo'])}–{r(q['crude_hi'])})",
           f"{r(q['mh_rr'])} ({r(q['mh_lo'])}–{r(q['mh_hi'])})"])

# ---------------------------------------------------------------- Table 3
t3 = [['', 'Admissions, n', 'Operative, n (%)', 'Operative mortality, %',
       'Nonoperative mortality, %']]
order = {'beds': ['<300', '300-599', '600-999', '1000+'],
         'owner': ['Oeffentlich', 'Freigemeinnuetzig', 'Privat']}
name = {'<300': '   <300', '300-599': '   300–599', '600-999': '   600–999',
        '1000+': '   ≥1000', 'Oeffentlich': '   Public',
        'Freigemeinnuetzig': '   Non-profit', 'Privat': '   Private'}
hx = {(x['dimension'], x['category']): x for x in hosp}
for dim, head in (('beds', 'Bed-size class'), ('owner', 'Ownership')):
    t3.append([head, '', '', '', ''])
    for cat in order[dim]:
        x = hx[(dim, cat)]
        t3.append([name[cat], thousands(x['N_total']),
                   f"{thousands(x['N_operative'])} ({r(x['operative_share_pct'], 1)})",
                   r(x['mort_operative_pct']), r(x['mort_nonoperative_pct'])])

# ------------------------------------------------- supplementary table content
s_age = [['Age band, years', 'Admissions, n', 'Operative, n (%)',
          'Operative mortality, %', 'Nonoperative mortality, %']]
for x in age:
    lab = x['label'].replace('-', '–').replace('80+', '≥80')
    s_age.append([lab, thousands(x['N_total']),
                  f"{thousands(x['N_operative'])} ({r(x['operative_share_pct'], 1)})",
                  r(x['mort_operative_pct']), r(x['mort_nonoperative_pct'])])

s_annual = [['Data year', 'Admissions, n', 'Operative, n (%)', 'Male, %',
             'Aged ≥80 years, %', 'Mean length of stay, days', 'PCCL 4–6, %',
             'Intensive care, %', 'Ventilation >96 h, %', 'Deaths, n (%)']]
for y in ann:
    a = ann[y]
    s_annual.append([y, thousands(a['episodes']),
                     f"{thousands(a['operative'])} ({r(a['operative_share_pct'], 1)})",
                     r(a['male_pct'], 1), r(a['age80plus_pct'], 1), r(a['mean_los_days'], 1),
                     r(a['pccl4_6_pct'], 1), r(a['icu_pct'], 1), r(a['vent96_pct'], 1),
                     f"{thousands(a['deaths'])} ({r(a['mortality_pct'])})"])

s_sens = [['Exposure definition', 'Admissions, n', 'Deaths, n', 'Mortality, %',
           'Crude RR versus nonoperative (95% CI)']]
for x in sens:
    s_sens.append([x['definition'].replace('core (primary)', 'Core (primary)')
                   .replace('wide', 'Wide').replace('instrumented only', 'Instrumented only'),
                   thousands(x['N_exposed']), thousands(x['deaths_exposed']),
                   r(x['mortality_pct']),
                   f"{r(x['rr'])} ({r(x['lo'])}–{r(x['hi'])})"])

fam_path = os.path.join(RES, 'procedure_families.csv')
s_fam = []
if os.path.exists(fam_path):
    s_fam = [['OPS family', 'Terminal codes per year', 'Episodes 2019-2024, n',
              '% of the operative arm']]
    for x in rd('procedure_families.csv'):
        codes = ('' if x['terminal_codes_min'] == ''
                 else (x['terminal_codes_min'] if x['terminal_codes_min'] == x['terminal_codes_max']
                       else x['terminal_codes_min'] + '-' + x['terminal_codes_max']))
        s_fam.append([x['family'], codes, thousands(x['episodes_2019_2024']),
                      r(x['pct_of_operative_arm'], 1)])

coh_path = os.path.join(RES, 'cohort_definition.csv')
s_coh = []
if os.path.exists(coh_path):
    s_coh = [['Cohort definition', 'Admissions, n', 'Operative, %', 'Mortality, %',
              'Operative mortality, %', 'Nonoperative mortality, %', 'Crude RR (95% CI)']]
    LAB = {'principal diagnosis (primary cohort)': 'Principal diagnosis (primary cohort)',
           'secondary diagnosis (overlap included)': 'Secondary diagnosis, overlap included',
           'secondary diagnosis only (overlap removed)': 'Secondary diagnosis only',
           'overlap (both principal and secondary)': 'Overlap (principal and secondary)',
           'any diagnosis (union)': 'Any diagnosis (union)',
           'secondary diagnosis (overlap included), age-adjusted':
               'Secondary diagnosis, overlap included, age-adjusted'}
    for x in rd('cohort_definition.csv'):
        s_coh.append([LAB.get(x['cohort'], x['cohort']), thousands(x['N_total']),
                      r(x['operative_share_pct'], 1), r(x['mortality_pct']),
                      r(x['mort_operative_pct']), r(x['mort_nonoperative_pct']),
                      f"{r(x['rr'])} ({r(x['lo'])}–{r(x['hi'])})"])

spec_path = os.path.join(RES, 'diagnosis_specificity.csv')
s_spec = []
if os.path.exists(spec_path):
    s_spec = [['Diagnosis code group', 'Admissions, n', 'Share of the cohort, %', 'Deaths, n',
               'Mortality, %', 'Mortality RR (95% CI)']]
    for x in rd('diagnosis_specificity.csv'):
        rrtxt = ('1 (reference)' if x['lo'] == ''
                 else f"{r(x['rr_vs_unspecified'])} ({r(x['lo'])}–{r(x['hi'])})")
        s_spec.append([x['code_group'], thousands(x['episodes']), r(x['share_of_cohort_pct'], 1),
                       thousands(x['deaths']), r(x['mortality_pct']), rrtxt])

ind_path = os.path.join(RES, 'indication_comorbidity.csv')
s_ind = []
if os.path.exists(ind_path):
    s_ind = [['Condition (ICD-10-GM secondary diagnosis)', 'Operative, n (%)',
              'Nonoperative, n (%)', 'Prevalence ratio (95% CI)']]
    for x in rd('indication_comorbidity.csv'):
        if x['prevalence_ratio'] == '':
            s_ind.append([x['condition'], r(x['pct_operative']), r(x['pct_nonoperative']), '—'])
            continue
        s_ind.append([x['condition'],
                      f"{thousands(x['N_operative'])} ({r(x['pct_operative'])})",
                      f"{thousands(x['N_nonoperative'])} ({r(x['pct_nonoperative'])})",
                      f"{r(x['prevalence_ratio'])} ({r(x['lo'])}–{r(x['hi'])})"])

adj_path = os.path.join(RES, 'adjustment_sets.csv')
s_adj = []
if os.path.exists(adj_path):
    s_adj = [['Stratification variables', 'Risk ratio (95% CI)', 'Strata']]
    for x in rd('adjustment_sets.csv'):
        lab = x['adjustment'][0].upper() + x['adjustment'][1:]
        s_adj.append([lab, f"{r(x['rr'])} ({r(x['lo'])}–{r(x['hi'])})", x['note']])

dis_path = os.path.join(RES, 'discharge_destination.csv')
s_dis = []
if os.path.exists(dis_path):
    s_dis = [['Discharge destination', 'Operative, n (%)', 'Nonoperative, n (%)',
              'Whole cohort, n (%)']]
    for x in rd('discharge_destination.csv'):
        if x['N_operative'] == '':
            s_dis.append([x['destination'], 'Suppressed', 'Suppressed',
                          f"{thousands(x['N_cohort'])} ({r(x['pct_cohort'])})"])
            continue
        s_dis.append([x['destination'],
                      f"{thousands(x['N_operative'])} ({r(x['pct_operative'])})",
                      f"{thousands(x['N_nonoperative'])} ({r(x['pct_nonoperative'])})",
                      f"{thousands(x['N_cohort'])} ({r(x['pct_cohort'])})"])

tra_path = os.path.join(RES, 'transfer_sensitivity.csv')
s_tra = []
if os.path.exists(tra_path):
    s_tra = [['Assumed mortality after transfer', 'Year-stratified risk ratio']]
    for x in rd('transfer_sensitivity.csv'):
        if 'nulls' in x['scenario']:
            s_tra.append([x['scenario'][0].upper() + x['scenario'][1:], r(x['rr'], 1) + '%'])
        else:
            s_tra.append([x['scenario'], r(x['rr'])])

bf_path = os.path.join(RES, 'bounding_factors.csv')
s_bf = []
if os.path.exists(bf_path):
    s_bf = [['Measured indication', 'Association with operative coding',
             'Association with in-hospital death', 'Bounding factor']]
    for x in rd('bounding_factors.csv'):
        s_bf.append([x['condition'][0].upper() + x['condition'][1:],
                     r(x['rr_exposure_confounder']), r(x['rr_confounder_outcome']),
                     r(x['bounding_factor'])])

s_smr = [['Data year', 'Observed deaths', 'Expected deaths', 'SMR (95% CI)']]
for y, x in smr.items():
    s_smr.append([y.capitalize() if y == 'pooled' else y, thousands(x['observed']),
                  r(x['expected'], 1), f"{r(x['smr'])} ({r(x['lo'])}–{r(x['hi'])})"])

out = {
    'table1': t1, 'table2': t2, 'table3': t3,
    's_age': s_age, 's_annual': s_annual, 's_sens': s_sens, 's_smr': s_smr, 's_fam': s_fam, 's_coh': s_coh, 's_spec': s_spec,
    's_ind': s_ind, 's_adj': s_adj, 's_dis': s_dis, 's_tra': s_tra, 's_bf': s_bf,
    'standardised': {
        'operative': (r(std['directly standardised mortality, operative, %']['estimate']),
                      r(std['directly standardised mortality, operative, %']['lo']),
                      r(std['directly standardised mortality, operative, %']['hi'])),
        'non_operative': (r(std['directly standardised mortality, non-operative, %']['estimate']),
                          r(std['directly standardised mortality, non-operative, %']['lo']),
                          r(std['directly standardised mortality, non-operative, %']['hi'])),
        'difference': (r(std['standardised difference, percentage points']['estimate']),
                       r(std['standardised difference, percentage points']['lo']),
                       r(std['standardised difference, percentage points']['hi'])),
        'Q': r(std['between-year heterogeneity of the age-adjusted RR, Cochran Q']['estimate']),
        'I2': r(std['between-year heterogeneity, I2, %']['estimate'], 0),
        'P': r(std['between-year heterogeneity, P value']['estimate']),
    },
    'envelope': rd('rounding_envelope.csv')[0],
}
for n, t in (('table1.csv', t1), ('table2.csv', t2), ('table3.csv', t3)):
    with open(os.path.join(RES, n), 'w', newline='') as f:
        csv.writer(f).writerows(t)
with open(os.path.join(RES, 'tables.json'), 'w') as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
print('tables written to', RES)
