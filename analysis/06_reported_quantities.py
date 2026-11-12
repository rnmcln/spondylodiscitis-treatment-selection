#!/usr/bin/env python3
"""Hospital-stratified risk ratios, Table 1 and 2 counts, figure data and key numbers.

Reads the derived results written by 01_analysis.py and 05b_hospital_reconstruction.py
(results/) and the per-year arm totals and distributions (data/arm_totals.csv,
data/arm_distributions.csv; SPONDYLO_DATA). Writes:

  results/hospital_structure_rr.csv      Mantel-Haenszel risk ratios stratified by data year
                                         and hospital bed-size category and, separately, by
                                         data year and ownership, from the reconstruction of
                                         the hidden hospital cells (05b)
  results/hospital_structure_detail.csv  per-category counts, operative share and mortality
                                         (reconstruction medians and ranges)
  results/table1_article.csv             arm characteristics as n (%), Table 1
  results/table2_article.csv             mortality comparison under each stratification, Table 2
  results/figure1_data.csv               operative share and mortality by age group, Figure 1
  results/figure2_data.csv               operative share and mortality by bed size, Figure 2
  results/key_numbers.json               the numbers quoted in the article text

Deterministic; no random component. Standard library only.
"""
import csv, json, math, os
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get('SPONDYLO_RESULTS', os.path.join(HERE, '..', 'results'))
DATA = os.environ.get('SPONDYLO_DATA', os.path.join(HERE, '..', 'data'))
Z = 1.959963984540054

rd = lambda p: list(csv.DictReader(open(p, encoding='utf-8')))


def write(name, rows):
    with open(os.path.join(RES, name), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def mh_rr(cells):
    """Mantel-Haenszel risk ratio with Greenland-Robins variance (as 01_analysis.py)."""
    num = den = P = 0.0
    for Ns, Ds, Nc, Dc in cells:
        T = Ns + Nc
        if T == 0 or Ns == 0 or Nc == 0:
            continue
        num += Ds * Nc / T
        den += Dc * Ns / T
        P += ((Ds + Dc) * Ns * Nc - Ds * Dc * T) / T**2
    rr = num / den
    se = math.sqrt(P / (num * den))
    return rr, rr * math.exp(-Z * se), rr * math.exp(Z * se)


def katz_rr(a, n1, b, n2):
    rr = (a / n1) / (b / n2)
    se = math.sqrt(1/a - 1/n1 + 1/b - 1/n2)
    return rr, rr * math.exp(-Z * se), rr * math.exp(Z * se)


def wilson(x, n):
    p = x / n
    d = 1 + Z**2 / n
    c = (p + Z**2 / (2*n)) / d
    h = Z * math.sqrt(p*(1-p)/n + Z**2/(4*n*n)) / d
    return 100*(c-h), 100*(c+h)


def pr_ci(a, n1, b, n2):
    return katz_rr(a, n1, b, n2)


# ------------------------------------------------------------ inputs
age = rd(os.path.join(RES, 'age_specific.csv'))
adj = {r['adjustment']: r for r in rd(os.path.join(RES, 'adjustment_sets.csv'))}
std = {r['quantity']: r for r in rd(os.path.join(RES, 'standardised_mortality.csv'))}
ind = rd(os.path.join(RES, 'indication_comorbidity.csv'))
dist = rd(os.path.join(DATA, 'arm_distributions.csv'))
cb = rd(os.path.join(DATA, 'arm_totals.csv'))
disc = {r['destination']: r for r in rd(os.path.join(RES, 'discharge_destination.csv'))}
annual = rd(os.path.join(RES, 'annual_summary.csv'))

# arm totals
N1 = sum(int(r['N_operative']) for r in cb)
N0 = sum(int(r['N_nonoperative']) for r in cb)
D1 = sum(int(r['deaths_operative']) for r in cb)
D0 = sum(int(r['deaths_nonoperative']) for r in cb)
assert (N1, N0, D1, D0) == (26130, 35004, 1841, 2207), (N1, N0, D1, D0)

# ------------------------------------------------ 1. hospital structure
# Taken from 05b_hospital_reconstruction.py, which reconstructs the hidden cells of
# the year x bed-size x ownership tables. Medians over 2000 feasible reconstructions
# are used as point values; the 2.5th-97.5th percentile and the full range describe
# the uncertainty due to the hidden cells.
REC = {r['quantity']: r for r in rd(os.path.join(RES, 'hospital_reconstruction_summary.csv'))}


def rec(q, col='median'):
    return float(REC[q][col])


hs_rows = []
for dim, lab in (('year_beds', 'data year and hospital bed-size category'),
                 ('year_owner', 'data year and hospital ownership')):
    hs_rows.append({'stratification': lab, 'strata': 24 if dim == 'year_beds' else 18,
                    'mh_rr': rec(f'rr_{dim}_ci_at_median_draw'), 'mh_lo': rec(f'rr_{dim}_ci_at_median_draw', 'p2_5'),
                    'mh_hi': rec(f'rr_{dim}_ci_at_median_draw', 'p97_5'),
                    'reconstruction_p2_5': rec(f'rr_{dim}', 'p2_5'), 'reconstruction_p97_5': rec(f'rr_{dim}', 'p97_5'),
                    'reconstruction_min': rec(f'rr_{dim}', 'min'), 'reconstruction_max': rec(f'rr_{dim}', 'max')})
write('hospital_structure_rr.csv', hs_rows)
hs = {'beds': hs_rows[0], 'owner': hs_rows[1]}

detail = []
for dim, cats in (('beds', ['<300', '300-599', '600-999', '1000+']),
                  ('owner', ['Oeffentlich', 'Freigemeinnuetzig', 'Privat'])):
    for c in cats:
        g = lambda f, col='median': rec(f'{dim}|{c}|{f}', col)
        n1, n2, dd, d1 = g('admissions'), g('operative'), g('deaths'), g('deaths_operative')
        detail.append({'dimension': dim, 'category': c, 'admissions': round(n1), 'operative': round(n2),
                       'operative_share_pct': round(g('operative_share_pct'), 1),
                       'share_range': f"{g('operative_share_pct', 'min'):.1f}-{g('operative_share_pct', 'max'):.1f}",
                       'deaths': round(dd), 'deaths_operative': round(d1),
                       'op_mort_pct': round(g('mort_operative_pct'), 2),
                       'nonop_mort_pct': round(g('mort_nonoperative_pct'), 2),
                       'stratum_rr': round(g('mort_operative_pct') / g('mort_nonoperative_pct'), 2)})
write('hospital_structure_detail.csv', detail)

# ------------------------------------------------ 2. Table 1 counts
def recon(state, col):
    """Sum over years of pct x N / 100 for a distribution column (reconstructed count)."""
    return sum(float(r[col]) * int(r['fallzahl']) / 100 for r in dist if r['state'] == state)

def recon_mean(state, col):
    return sum(float(r[col]) * int(r['fallzahl']) for r in dist if r['state'] == state)

male1 = recon('operative', 'male_pct'); male0 = recon('cohort', 'male_pct') - male1
pccl1 = sum(recon('operative', c) for c in ('pccl4', 'pccl5', 'pccl6'))
pccl0 = sum(recon('cohort', c) for c in ('pccl4', 'pccl5', 'pccl6')) - pccl1
los1 = recon_mean('operative', 'vwd_mean') / N1
los0 = (recon_mean('cohort', 'vwd_mean') - recon_mean('operative', 'vwd_mean')) / N0
icu1 = sum(int(r['ICU_operative']) for r in cb); icu0 = sum(int(r['ICU_total']) for r in cb) - icu1
ven1 = sum(int(r['vent96_operative']) for r in cb); ven0 = sum(int(r['vent96_total']) for r in cb) - ven1

ag = {r['band']: r for r in age}
def agesum(keys, col):
    return sum(int(ag[k][col]) for k in keys)
young = ['age_u30', 'age3039', 'age4049', 'age5054', 'age5559', 'age6064']
mid = ['age6574', 'age7579']
old = ['age80p']

ic = {r['condition'].split(' (')[0]: r for r in ind if r['prevalence_ratio']}
def cond(label, key):
    r = ic[key]
    a, b = int(r['N_operative']), int(r['N_nonoperative'])
    p = pr_ci(a, N1, b, N0)
    return (label, a, b, p)

t1 = []
def addrow(label, a, b, pr=None, kind='n'):
    row = OrderedDict(characteristic=label, kind=kind)
    if kind == 'n':
        row.update(op_n=round(a), op_pct=round(100*a/N1, 1), nonop_n=round(b), nonop_pct=round(100*b/N0, 1))
    else:
        row.update(op_n='', op_pct=round(a, 1), nonop_n='', nonop_pct=round(b, 1))
    if pr:
        row.update(pr=round(pr[0], 2), pr_lo=round(pr[1], 2), pr_hi=round(pr[2], 2))
    else:
        row.update(pr='', pr_lo='', pr_hi='')
    t1.append(row)

addrow('Male sex', male1, male0)
addrow('Age <65 years', agesum(young, 'N_operative'), agesum(young, 'N_nonoperative'))
addrow('Age 65-79 years', agesum(mid, 'N_operative'), agesum(mid, 'N_nonoperative'))
addrow('Age >=80 years', agesum(old, 'N_operative'), agesum(old, 'N_nonoperative'))
for label, key in [('Intraspinal empyema (G06.1)', 'Intraspinal / epidural empyema'),
                   ('Septic shock (R57.2)', 'Septic shock'),
                   ('Sepsis (A41.-)', 'Sepsis'),
                   ('Paraparesis or paraplegia (G82.-)', 'Paraparesis / paraplegia'),
                   ('Type 2 diabetes (E11.-)', 'Type 2 diabetes'),
                   ('Heart failure (I50.-)', 'Heart failure'),
                   ('Chronic kidney disease (N18.-)', 'Chronic kidney disease')]:
    l, a, b, p = cond(label, key)
    addrow(l, a, b, p)
addrow('PCCL 4-6', pccl1, pccl0)
addrow('Intensive care', icu1, icu0)
addrow('Mechanical ventilation >96 h', ven1, ven0)
addrow('Mean length of stay, days', los1, los0, kind='mean')
addrow('Transfer to another acute hospital',
       int(disc['Transfer to another acute hospital']['N_operative']),
       int(disc['Transfer to another acute hospital']['N_nonoperative']))
addrow('In-hospital death', D1, D0)
write('table1_article.csv', t1)

# ------------------------------------------------ 3. Table 2 rows
crude = katz_rr(D1, N1, D0, N0)
t2 = [
    ('None (crude)', 1, crude),
    ('Data year', 6, tuple(float(adj['data year'][k]) for k in ('rr', 'lo', 'hi'))),
    ('Data year and age group (primary analysis)', 54,
     tuple(float(adj['data year and age band'][k]) for k in ('rr', 'lo', 'hi'))),
    ('Data year and intraspinal empyema', 12,
     tuple(float(adj['data year and intraspinal empyema (G06.1)'][k]) for k in ('rr', 'lo', 'hi'))),
    ('Data year and sepsis', 12,
     tuple(float(adj['data year and sepsis (A41.-)'][k]) for k in ('rr', 'lo', 'hi'))),
    ('Data year and hospital bed-size category', 24, (hs['beds']['mh_rr'], hs['beds']['mh_lo'], hs['beds']['mh_hi'])),
    ('Data year and hospital ownership', 18, (hs['owner']['mh_rr'], hs['owner']['mh_lo'], hs['owner']['mh_hi'])),
]
write('table2_article.csv', [{'stratification': s, 'strata': k, 'rr': round(v[0], 4),
                          'lo': round(v[1], 4), 'hi': round(v[2], 4)} for s, k, v in t2])

# ------------------------------------------------ 4. figure data
f1 = []
for r in age:
    n1, d1 = int(r['N_operative']), int(r['D_operative'])
    n0, d0 = int(r['N_nonoperative']), int(r['D_nonoperative'])
    lo1, hi1 = wilson(d1, n1); lo0, hi0 = wilson(d0, n0)
    f1.append({'age_group': r['label'], 'admissions': int(r['N_total']), 'operative_n': n1,
               'operative_share_pct': float(r['operative_share_pct']),
               'op_deaths': d1, 'op_mort_pct': round(100*d1/n1, 2), 'op_lo': round(lo1, 2), 'op_hi': round(hi1, 2),
               'nonop_n': n0, 'nonop_deaths': d0, 'nonop_mort_pct': round(100*d0/n0, 2),
               'nonop_lo': round(lo0, 2), 'nonop_hi': round(hi0, 2)})
write('figure1_data.csv', f1)

order = ['<300', '300-599', '600-999', '1000+']
f2 = []
for k in order:
    g = lambda f, col='median': rec(f'beds|{k}|{f}', col)
    N, n1, D, d1 = g('admissions'), g('operative'), g('deaths'), g('deaths_operative')
    n0, d0 = N - n1, D - d1
    lo1, hi1 = wilson(d1, n1); lo0, hi0 = wilson(d0, n0); slo, shi = wilson(n1, N)
    f2.append({'bed_size': k, 'admissions': round(N), 'operative': round(n1),
               'operative_share_pct': round(g('operative_share_pct'), 2), 'share_lo': round(slo, 2),
               'share_hi': round(shi, 2), 'share_min': round(g('operative_share_pct', 'min'), 2),
               'share_max': round(g('operative_share_pct', 'max'), 2),
               'op_mort_pct': round(g('mort_operative_pct'), 2), 'op_lo': round(lo1, 2), 'op_hi': round(hi1, 2),
               'nonop_mort_pct': round(g('mort_nonoperative_pct'), 2), 'nonop_lo': round(lo0, 2),
               'nonop_hi': round(hi0, 2)})
write('figure2_data.csv', f2)

# ------------------------------------------------ 5. key numbers
key = OrderedDict()
key['admissions'] = N1 + N0
key['operative_n'] = N1
key['operative_pct'] = round(100*N1/(N1+N0), 1)
key['deaths'] = D1 + D0
key['mortality_pct'] = round(100*(D1+D0)/(N1+N0), 2)
key['op_mort_pct'] = round(100*D1/N1, 2)
key['nonop_mort_pct'] = round(100*D0/N0, 2)
key['crude_rr'] = [round(x, 2) for x in crude]
key['annual_admissions'] = {r['year']: int(r['episodes']) for r in annual}
key['grid_admissions'] = rec('grid_admissions')
key['outside_grid_admissions'] = rec('outside_grid_Q1')
key['share_admissions_lt300_pct'] = round(rec('share_admissions_lt300_pct'), 1)
key['share_operations_ge600_pct'] = round(rec('share_operations_ge600_pct'), 1)
key['op_share_by_bed'] = {k: round(rec(f'beds|{k}|operative_share_pct'), 1) for k in order}
key['op_share_by_age'] = {r['label']: float(r['operative_share_pct']) for r in age}
key['beds_rr'] = [hs['beds'][k] for k in ('mh_rr', 'mh_lo', 'mh_hi')]
key['beds_rr_range'] = [hs['beds']['reconstruction_min'], hs['beds']['reconstruction_max']]
key['owner_rr_range'] = [hs['owner']['reconstruction_min'], hs['owner']['reconstruction_max']]
key['owner_rr'] = [hs['owner'][k] for k in ('mh_rr', 'mh_lo', 'mh_hi')]
key['age_adjusted_rr'] = [float(adj['data year and age band'][k]) for k in ('rr', 'lo', 'hi')]
key['std_op'] = float(std['directly standardised mortality, operative, %']['estimate'])
key['std_nonop'] = float(std['directly standardised mortality, non-operative, %']['estimate'])
key['std_diff'] = [float(std['standardised difference, percentage points'][k]) for k in ('estimate', 'lo', 'hi')]
key['table1'] = t1
key['figure2'] = f2
with open(os.path.join(RES, 'key_numbers.json'), 'w', encoding='utf-8') as fh:
    json.dump(key, fh, indent=1, ensure_ascii=False)

for k, v in key.items():
    if k not in ('table1', 'figure2'):
        print(f'{k}: {v}')
print()
for r in t1:
    print(r)
print()
for r in hs_rows:
    print(r)
