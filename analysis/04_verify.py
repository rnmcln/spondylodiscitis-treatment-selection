#!/usr/bin/env python3
"""Consistency checks on the input data and the derived results.

Run after 01_analysis.py. Exits non-zero if any check fails. The checks cover
transcription integrity of the query output, the internal consistency of the
derived arms, and agreement between the results files and the headline values
reported in the paper.
"""
import csv, math, os, sys
from decimal import Decimal, ROUND_HALF_UP

DATA = os.environ.get('SPONDYLO_DATA', os.path.join(os.path.dirname(__file__), '..', 'data'))
RES = os.environ.get('SPONDYLO_RESULTS', os.path.join(os.path.dirname(__file__), '..', 'results'))
AGE = ['age_u30', 'age3039', 'age4049', 'age5054', 'age5559', 'age6064',
       'age6574', 'age7579', 'age80p']
PCCL = ['pccl0', 'pccl1', 'pccl2', 'pccl3', 'pccl4', 'pccl5', 'pccl6']

rdd = lambda n: list(csv.DictReader(open(os.path.join(DATA, n))))
rdr = lambda n: list(csv.DictReader(open(os.path.join(RES, n))))
r2 = lambda x: str(Decimal(repr(float(x))).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP))

fails = []
n = 0


def check(label, ok, detail=''):
    global n
    n += 1
    print(('PASS ' if ok else 'FAIL ') + label + ((' [' + detail + ']') if detail and not ok else ''))
    if not ok:
        fails.append(label)


totals = {int(r['year']): r for r in rdd('arm_totals.csv')}
dist = {(int(r['year']), r['state']): r for r in rdd('arm_distributions.csv')}
kv = rdd('hospital_tables.csv')
years = sorted(totals)

print('--- input integrity ---')
for y in years:
    t = totals[y]
    check(f'{y}: operative + non-operative = total',
          int(t['N_operative']) + int(t['N_nonoperative']) == int(t['N_total']))
    check(f'{y}: operative + non-operative deaths = total deaths',
          int(t['deaths_operative']) + int(t['deaths_nonoperative']) == int(t['deaths_total']))
    check(f'{y}: deaths do not exceed episodes in either arm',
          int(t['deaths_operative']) <= int(t['N_operative']) and
          int(t['deaths_nonoperative']) <= int(t['N_nonoperative']))
    check(f'{y}: instrumented subset does not exceed the core exposure',
          int(t['N_instrumented']) <= int(t['N_operative']))
    check(f'{y}: wide definition is not smaller than the core exposure',
          int(t['N_operative_wide']) >= int(t['N_operative']))
    for state, key in (('cohort', 'N_total'), ('operative', 'N_operative'),
                       ('deaths_total', 'deaths_total'), ('deaths_operative', 'deaths_operative')):
        d = dist[(y, state)]
        check(f'{y} {state}: distribution row matches the query total',
              int(float(d['fallzahl'])) == int(t[key]))
        check(f'{y} {state}: age distribution sums to 100 within rounding',
              abs(sum(float(d[k]) for k in AGE) - 100) < 0.09)
        check(f'{y} {state}: PCCL distribution sums to 100',
              abs(sum(float(d[k]) for k in PCCL) - 100) < 0.05)

for y in years:
    for q, key in (('Q1', 'N_total'), ('Q2', 'N_operative'),
                   ('Q3', 'deaths_total'), ('Q4', 'deaths_operative')):
        s = sum(float(r['faelle']) for r in kv
                if int(r['year']) == y and r['query'] == q)
        check(f'{y} {q}: hospital table sums to the query total',
              abs(s - int(totals[y][key])) < 0.5)

print('\n--- derived results ---')
rr = {r['year']: r for r in rdr('risk_ratios.csv')}
age = {r['band']: r for r in rdr('age_specific.csv')}
std = {r['quantity']: r for r in rdr('standardised_mortality.csv')}
sens = {r['definition']: r for r in rdr('exposure_sensitivity.csv')}
env = rdr('rounding_envelope.csv')[0]

g = lambda k: sum(int(totals[y][k]) for y in years)
N, O = g('N_total'), g('N_operative')
D, DO = g('deaths_total'), g('deaths_operative')
check('pooled episodes = 61,134', N == 61134)
check('pooled operative = 26,130', O == 26130)
check('pooled deaths = 4048', D == 4048)
check('operative mortality = 7.05%', r2(100*DO/O) == '7.05')
check('non-operative mortality = 6.30%', r2(100*(D-DO)/(N-O)) == '6.30')
check('pooled crude RR = 1.12 (1.05-1.19)',
      (r2(rr['pooled']['crude_rr']), r2(rr['pooled']['crude_lo']), r2(rr['pooled']['crude_hi']))
      == ('1.12', '1.05', '1.19'))
check('pooled age-adjusted RR = 1.31 (1.24-1.39)',
      (r2(rr['pooled']['mh_rr']), r2(rr['pooled']['mh_lo']), r2(rr['pooled']['mh_hi']))
      == ('1.31', '1.24', '1.39'))
check('standardised difference = 1.84 pp (1.43-2.25)',
      (r2(std['standardised difference, percentage points']['estimate']),
       r2(std['standardised difference, percentage points']['lo']),
       r2(std['standardised difference, percentage points']['hi'])) == ('1.84', '1.43', '2.25'))
check('operative share falls monotonically from age 50-54 upwards',
      all(float(age[a]['operative_share_pct']) > float(age[b]['operative_share_pct'])
          for a, b in zip(['age5054', 'age5559', 'age6064', 'age6574', 'age7579'],
                          ['age5559', 'age6064', 'age6574', 'age7579', 'age80p'])))
check('operative mortality exceeds non-operative in every band above 30',
      all(float(age[b]['mort_operative_pct']) > float(age[b]['mort_nonoperative_pct'])
          for b in AGE[1:]))
check('every exposure definition gives a risk ratio above 1',
      all(float(v['rr']) > 1 and float(v['lo']) > 1 for v in sens.values()))
check('rounding envelope spans less than 0.01 of the risk ratio',
      float(env['max']) - float(env['min']) < 0.01)
off = {(y, st) for (y, st), d in dist.items()
       if abs(sum(float(d[k]) for k in AGE) - 100) > 0.09}
check('every age distribution sums to 100 within the rounding of sixteen published bands',
      not off, str(sorted(off)))
check('E-value of the pooled estimate = 1.95',
      r2(float(rr['pooled']['mh_rr']) +
         math.sqrt(float(rr['pooled']['mh_rr'])*(float(rr['pooled']['mh_rr'])-1))) == '1.95')


fam_path = os.path.join(RES, 'procedure_families.csv')
if os.path.exists(fam_path):
    fam = {r['family']: r for r in rdr('procedure_families.csv')}
    check('every operative family is smaller than the operative arm',
          all(int(v['episodes_2019_2024']) < O for k, v in fam.items() if k != 'sum over families'))
    check('families overlap: their sum exceeds the operative arm',
          int(fam['sum over families']['episodes_2019_2024']) > O)
    fc_path = os.path.join(DATA, 'procedure_family_counts.csv')
    if os.path.exists(fc_path):
        rows = rdd('procedure_family_counts.csv')
        by_year = {}
        for r in rows:
            by_year[r['year']] = by_year.get(r['year'], 0) + int(r['codes'])
        expected = {'2019': 174, '2020': 180, '2021': 218, '2022': 218, '2023': 212, '2024': 214}
        check('per-family terminal codes sum to the core code count in every year',
              by_year == expected, str(by_year))

coh_path = os.path.join(RES, 'cohort_definition.csv')
if os.path.exists(coh_path):
    coh = {r['cohort']: r for r in rdr('cohort_definition.csv')}
    check('the secondary-diagnosis cohort reverses the crude association',
          float(coh['secondary diagnosis (overlap included)']['rr']) < 1 <
          float(coh['principal diagnosis (primary cohort)']['rr']))
    uni = coh.get('any diagnosis (union)')
    ovl = coh.get('overlap (both principal and secondary)')
    if uni and ovl:
        check('the union equals principal plus secondary minus the measured overlap',
              int(uni['N_total']) == int(coh['principal diagnosis (primary cohort)']['N_total'])
              + int(coh['secondary diagnosis (overlap included)']['N_total']) - int(ovl['N_total']))
        check('the overlap is a subset of both cohorts',
              int(ovl['N_total']) < int(coh['principal diagnosis (primary cohort)']['N_total'])
              and int(ovl['N_total']) < int(coh['secondary diagnosis (overlap included)']['N_total']))
    else:
        check('the pooled any-diagnosis cohort is the sum of the two',
              int(coh['any diagnosis (sum of the two, overlap not measured)']['N_total']) ==
              int(coh['principal diagnosis (primary cohort)']['N_total']) +
              int(coh['secondary diagnosis (overlap included)']['N_total']))
    adj = coh.get('secondary diagnosis (overlap included), age-adjusted')
    if adj:
        check('age adjustment moves the secondary-diagnosis estimate towards the null',
              float(coh['secondary diagnosis (overlap included)']['rr']) < float(adj['rr']) < 1)
        check('age adjustment moves the primary-cohort estimate away from the null',
              float(rr['pooled']['crude_rr']) < float(rr['pooled']['mh_rr']))
        sd_path = os.path.join(DATA, 'secondary_diagnosis_distributions.csv')
        if os.path.exists(sd_path):
            sd = rdd('secondary_diagnosis_distributions.csv')
            check('secondary-diagnosis age distributions sum to 100 within rounding',
                  all(abs(sum(float(r[k]) for k in AGE) - 100) < 0.09 for r in sd))
            check('secondary-diagnosis distribution rows match their query totals',
                  len(sd) == 24)

coh_path = os.path.join(RES, 'cohort_definition.csv')
if os.path.exists(coh_path):
    cohr = {r['cohort']: r for r in rdr('cohort_definition.csv')}
    if 'secondary diagnosis only (overlap removed)' in cohr:
        only = cohr['secondary diagnosis only (overlap removed)']
        sec_ = cohr['secondary diagnosis (overlap included)']
        ovl_ = cohr['overlap (both principal and secondary)']
        uni_ = cohr['any diagnosis (union)']
        pri_ = cohr['principal diagnosis (primary cohort)']
        check('the secondary-only cohort is the secondary cohort minus the overlap',
              int(only['N_total']) == int(sec_['N_total']) - int(ovl_['N_total'])
              and int(only['N_operative']) == int(sec_['N_operative']) - int(ovl_['N_operative']))
        check('the union is the principal cohort plus the secondary-only cohort',
              int(uni_['N_total']) == int(pri_['N_total']) + int(only['N_total']))
        check('every cohort has a non-negative non-operative arm',
              all(int(r['N_total']) >= int(r['N_operative']) for r in cohr.values()))

db_path = os.path.join(DATA, 'exposure_debridement.csv')
if os.path.exists(db_path):
    dbr = {int(r['year']): r for r in rdd('exposure_debridement.csv')}
    check('the debridement-directed series covers every data year', sorted(dbr) == years)
    check('debridement-directed counts do not exceed the core operative arm',
          all(int(dbr[y]['N_debridement']) <= int(totals[y]['N_operative'])
              and int(dbr[y]['deaths_debridement']) <= int(totals[y]['deaths_operative'])
              for y in years))
    check('the debridement-directed definition selects fewer codes than the core definition',
          all(int(dbr[y]['codes']) < 174 for y in years))

spec_path = os.path.join(RES, 'diagnosis_specificity.csv')
if os.path.exists(spec_path):
    sp = {r['code_group']: r for r in rdr('diagnosis_specificity.csv')}
    a = sp['specific (M46.2, M46.3)']; b = sp['unspecified (M46.4, M46.5)']
    check('the diagnosis-specificity split reconciles with the cohort',
          int(a['episodes']) + int(b['episodes']) == N and
          int(a['deaths']) + int(b['deaths']) == D)
    check('the unspecified codes dominate the cohort',
          float(b['share_of_cohort_pct']) > 90)
    ds_path = os.path.join(DATA, 'diagnosis_specificity.csv')
    if os.path.exists(ds_path):
        rows = rdd('diagnosis_specificity.csv')
        check('every year of the specificity split sums to its cohort total',
              all(int(r['N_specific']) + int(r['N_unspecific']) ==
                  int(totals[int(r['year'])]['N_total']) for r in rows))
        check('every year of the specificity split sums to its death total',
              all(int(r['deaths_specific']) + int(r['deaths_unspecific']) ==
                  int(totals[int(r['year'])]['deaths_total']) for r in rows))

panel_path = os.path.join(DATA, 'secondary_diagnosis_panel.csv')
if os.path.exists(panel_path):
    print('--- secondary-diagnosis register, comorbidity and discharge inputs ---')
    pan = {(int(r['year']), r['q']): r for r in rdd('secondary_diagnosis_panel.csv')}
    check('the register was read for all four query states in every data year',
          len(pan) == 24)
    check('each register read reports the case count of its own query',
          all(int(pan[(y, q)]['n']) in (int(totals[y][k]) - 1, int(totals[y][k]),
                                        int(totals[y][k]) + 1)
              for y in years
              for q, k in (('Q1', 'N_total'), ('Q2', 'N_operative'),
                           ('Q3', 'deaths_total'), ('Q4', 'deaths_operative'))),
          'implied N is reconstructed from a percentage printed to two decimals')
    codes = [c for c in pan[(years[0], 'Q1')] if c not in ('year', 'q', 'n', 'nenn',
                                                           'meanDx', 'restF')]
    check('every register count in the operative arm is at most the cohort count',
          all(float(pan[(y, 'Q2')][c]) <= float(pan[(y, 'Q1')][c])
              for y in years for c in codes))
    check('every register count among operative deaths is at most that among all deaths',
          all(float(pan[(y, 'Q4')][c]) <= float(pan[(y, 'Q3')][c])
              for y in years for c in codes))
    check('coded secondary diagnosis entries per episode exceed ten in every query',
          all(float(pan[(y, q)]['meanDx']) > 10 for y in years
              for q in ('Q1', 'Q2', 'Q3', 'Q4')))
    if 'faelle' in pan[(years[0], 'Q1')]:
        check('entries are at least as numerous as distinct codes in every query',
              all(int(pan[(y, q)]['nenn']) >= int(pan[(y, q)]['faelle'])
                  for y in years for q in ('Q1', 'Q2', 'Q3', 'Q4')))
        check('the register lists between 300 and 2000 codes per query',
              all(300 <= int(pan[(y, q)]['codes']) <= 2000
                  for y in years for q in ('Q1', 'Q2', 'Q3', 'Q4')))
        check('the suppressed residual is counted once in the distinct-code sum',
              all(int(pan[(y, q)]['restFaelle']) == int(pan[(y, q)]['restF'])
                  for y in years for q in ('Q1', 'Q2', 'Q3', 'Q4')))

ex_path = os.path.join(DATA, 'comorbidity_exact.csv')
if os.path.exists(ex_path):
    exr = {(int(r['year']), r['arm'], r['group'], r['outcome']): int(r['n'])
           for r in rdd('comorbidity_exact.csv')}
    groups = sorted({k[2] for k in exr})
    check('the comorbidity grid is complete',
          len(exr) == len(years) * 2 * (len(groups) + 1))
    check('every operative count is at most the cohort count',
          all(exr[(y, 'operative', g, 'all')] <= exr[(y, 'cohort', g, 'all')]
              for y in years for g in groups))
    check('every comorbidity count is at most its arm total',
          all(exr[(y, 'operative', g, 'all')] <= int(totals[y]['N_operative'])
              and exr[(y, 'cohort', g, 'all')] <= int(totals[y]['N_total'])
              for y in years for g in groups))
    check('sepsis deaths are at most sepsis episodes and at most all deaths',
          all(exr[(y, a, 'A41', 'death')] <= exr[(y, a, 'A41', 'all')]
              for y in years for a in ('cohort', 'operative'))
          and all(exr[(y, 'cohort', 'A41', 'death')] <= int(totals[y]['deaths_total'])
                  and exr[(y, 'operative', 'A41', 'death')] <=
                  int(totals[y]['deaths_operative']) for y in years))

dd_path = os.path.join(DATA, 'discharge_destination.csv')
if os.path.exists(dd_path):
    ddr = {(int(r['year']), r['arm'], r['destination']): r['n']
           for r in rdd('discharge_destination.csv')}
    dests = sorted({k[2] for k in ddr})
    check('the discharge grid is complete', len(ddr) == len(years) * 2 * len(dests))
    check('every operative destination count is at most the cohort count',
          all(ddr[(y, 'operative', d)] == 'SUP' or ddr[(y, 'cohort', d)] == 'SUP'
              or int(ddr[(y, 'operative', d)]) <= int(ddr[(y, 'cohort', d)])
              for y in years for d in dests))
    check('destination counts and deaths do not exceed the arm totals',
          all(sum(int(ddr[(y, 'cohort', d)]) for d in dests
                  if ddr[(y, 'cohort', d)] != 'SUP')
              + int(totals[y]['deaths_total']) <= int(totals[y]['N_total'])
              for y in years))

print(f'\n{n} checks, {len(fails)} failures')
if fails:
    print('FAILED:', fails)
    sys.exit(1)
