#!/usr/bin/env python3
"""Treatment selection and in-hospital mortality in pyogenic spondylodiscitis.

Analysis of InEK DatenBrowser aggregate query results (principal diagnosis
M46.2-M46.5, data years 2019-2024).

Inputs (see analysis/schema/input_schema.md; regenerate them from the public
DatenBrowser using protocol/query_grid.csv):
  arm_totals.csv          per-year arm totals, sensitivity-definition totals,
                          complexity and resource-use marker counts
  arm_distributions.csv   per-year age, sex, PCCL and length-of-stay
                          distributions for each query state
  hospital_tables.csv     bed-size x ownership tables for Q1-Q4

Outputs are written to results/.

Analytical steps
  1. Non-operative quantities by subtraction (Q1-Q2, Q3-Q4).
  2. Age-stratum counts reconstructed from the published percentage
     distributions (two decimal places) and rescaled proportionally to the
     exact query totals.
  3. Operative share by age band and by hospital characteristics.
  4. Crude risk ratios (Katz logarithmic intervals); Mantel-Haenszel risk
     ratios over age strata with Greenland-Robins variance.
  5. Direct age standardisation of arm-specific mortality to the pooled cohort
     age distribution, with the standardised absolute difference.
  6. Between-year heterogeneity of the age-adjusted association (Cochran Q, I2).
  7. Rounding-envelope sensitivity: percentages resampled within their
     two-decimal rounding intervals, Mantel-Haenszel estimate recomputed.
  8. Exposure-definition sensitivity analyses; indirect standardisation (SMR)
     and E-values, reported as supplementary quantities.
"""
import csv, math, os, random
from collections import defaultdict

DATA = os.environ.get('SPONDYLO_DATA', os.path.join(os.path.dirname(__file__), '..', 'data'))
OUT = os.environ.get('SPONDYLO_RESULTS', os.path.join(os.path.dirname(__file__), '..', 'results'))
os.makedirs(OUT, exist_ok=True)

AGE_KEYS = ['age_u30', 'age3039', 'age4049', 'age5054', 'age5559', 'age6064',
            'age6574', 'age7579', 'age80p']
AGE_LABEL = {'age_u30': '<30', 'age3039': '30-39', 'age4049': '40-49', 'age5054': '50-54',
             'age5559': '55-59', 'age6064': '60-64', 'age6574': '65-74',
             'age7579': '75-79', 'age80p': '80+'}
Z = 1.959963984540054

rd = lambda n: list(csv.DictReader(open(os.path.join(DATA, n))))
res = {int(r['year']): r for r in rd('arm_totals.csv')}
dist = {(int(r['year']), r['state']): r for r in rd('arm_distributions.csv')}
kv = rd('hospital_tables.csv')
years = sorted(res)

# Optional inputs for the coded-condition analyses. Loaded here because Table 1
# reports three of the conditions alongside the other arm characteristics.
_pp = os.path.join(DATA, 'secondary_diagnosis_panel.csv')
_ep = os.path.join(DATA, 'comorbidity_exact.csv')
panel_in = {(int(r['year']), r['q']): r for r in rd('secondary_diagnosis_panel.csv')} \
    if os.path.exists(_pp) else None
exact_in = {(int(r['year']), r['arm'], r['group'], r['outcome']): int(r['n'])
            for r in rd('comorbidity_exact.csv')} if os.path.exists(_ep) else None


def age_counts(year, state, total):
    """Reconstruct stratum counts from percentages, rescaled to the exact total."""
    d = dist[(year, state)]
    raw = [float(d[k]) / 100.0 * total for k in AGE_KEYS]
    s = sum(raw)
    return [x * total / s for x in raw] if s > 0 else raw


def build_strata(perturb=None):
    """Return [(year, band, N_op, D_op, N_nonop, D_nonop)] and the clipping log.

    perturb: optional function(percentage) -> percentage, used by the
    rounding-envelope sensitivity analysis.
    """
    strata, clipped = [], []
    for y in years:
        r = res[y]
        Nt, Ns = int(r['N_total']), int(r['N_operative'])
        Dt, Ds = int(r['deaths_total']), int(r['deaths_operative'])
        if perturb is None:
            ct = age_counts(y, 'cohort', Nt);        cs = age_counts(y, 'operative', Ns)
            dt = age_counts(y, 'deaths_total', Dt);  ds = age_counts(y, 'deaths_operative', Ds)
        else:
            def pc(state, total):
                d = dist[(y, state)]
                raw = [max(0.0, perturb(float(d[k]))) / 100.0 * total for k in AGE_KEYS]
                s = sum(raw)
                return [x * total / s for x in raw] if s > 0 else raw
            ct = pc('cohort', Nt);       cs = pc('operative', Ns)
            dt = pc('deaths_total', Dt); ds = pc('deaths_operative', Ds)
        for i, band in enumerate(AGE_KEYS):
            Ns_a, Ds_a = cs[i], ds[i]
            Nc_a, Dc_a = ct[i] - cs[i], dt[i] - ds[i]
            if Nc_a < 0 or Dc_a < 0:
                clipped.append((y, band, round(Nc_a, 2), round(Dc_a, 2)))
                Nc_a, Dc_a = max(Nc_a, 0.0), max(Dc_a, 0.0)
            Ds_a = min(Ds_a, Ns_a); Dc_a = min(Dc_a, Nc_a)
            strata.append((y, band, Ns_a, Ds_a, Nc_a, Dc_a))
    return strata, clipped


def mh_rr(cells):
    """Mantel-Haenszel risk ratio with Greenland-Robins variance."""
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
    return rr, rr * math.exp(-Z * se), rr * math.exp(Z * se), se


def katz_rr(a, n1, b, n2):
    """Crude risk ratio with Katz logarithmic interval."""
    rr = (a / n1) / (b / n2)
    se = math.sqrt(1/a - 1/n1 + 1/b - 1/n2)
    return rr, rr * math.exp(-Z * se), rr * math.exp(Z * se)


def chi2_sf(x, k):
    """Upper tail probability of a chi-square distribution, standard library only.

    Regularised upper incomplete gamma Q(k/2, x/2), evaluated by the series
    expansion below the transition point and by the continued fraction above it.
    """
    a, xx = k/2.0, x/2.0
    if xx <= 0:
        return 1.0
    lg = math.lgamma(a)
    if xx < a + 1:                       # series for P(a, x), then Q = 1 - P
        term = 1.0/a
        s = term
        n = a
        for _ in range(1000):
            n += 1
            term *= xx/n
            s += term
            if abs(term) < abs(s)*1e-15:
                break
        return 1.0 - s*math.exp(-xx + a*math.log(xx) - lg)
    b, c = xx + 1 - a, 1e308             # continued fraction for Q(a, x)
    dd, h = 1.0/b, 1.0/b
    for i in range(1, 1000):
        an = -i*(i - a)
        b += 2
        dd = an*dd + b
        if abs(dd) < 1e-300:
            dd = 1e-300
        c = b + an/c
        if abs(c) < 1e-300:
            c = 1e-300
        dd = 1.0/dd
        delta = dd*c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return h*math.exp(-xx + a*math.log(xx) - lg)


def e_value(rr):
    if rr < 1:
        rr = 1 / rr
    return rr + math.sqrt(rr * (rr - 1))


def write(name, rows):
    with open(os.path.join(OUT, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f'  -> results/{name}')


strata, clipped = build_strata()
g = lambda k: sum(int(res[y][k]) for y in years)
Nt, Ns = g('N_total'), g('N_operative')
Dt, Ds = g('deaths_total'), g('deaths_operative')
Nc, Dc = Nt - Ns, Dt - Ds

print('=== Cohort and treatment distribution ===')
print(f'episodes {Nt}; operative {Ns} ({100*Ns/Nt:.1f}%); non-operative {Nc}')
print(f'deaths {Dt} ({100*Dt/Nt:.2f}%); operative {Ds} ({100*Ds/Ns:.2f}%); '
      f'non-operative {Dc} ({100*Dc/Nc:.2f}%)')

# ---------------------------------------------------------------- 1. arm profile
pool = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
for y, band, Nso, Dso, Nco, Dco in strata:
    p = pool[band]
    p[0] += Nso; p[1] += Dso; p[2] += Nco; p[3] += Dco

print('\n=== Operative share, mortality and age (pooled 2019-2024) ===')
age_rows = []
for band in AGE_KEYS:
    Nso, Dso, Nco, Dco = pool[band]
    tot = Nso + Nco
    ms = 100*Dso/Nso if Nso else float('nan')
    mc = 100*Dco/Nco if Nco else float('nan')
    print(f'{AGE_LABEL[band]:>6}: operative share {100*Nso/tot:5.1f}%   '
          f'mortality operative {ms:5.2f}%   non-operative {mc:5.2f}%')
    age_rows.append({'band': band, 'label': AGE_LABEL[band], 'N_total': round(tot),
                     'N_operative': round(Nso), 'operative_share_pct': round(100*Nso/tot, 1),
                     'D_operative': round(Dso), 'mort_operative_pct': round(ms, 2),
                     'N_nonoperative': round(Nco), 'D_nonoperative': round(Dco),
                     'mort_nonoperative_pct': round(mc, 2)})
write('age_specific.csv', age_rows)

def wsum(state, key):
    return sum(float(dist[(y, state)]['fallzahl']) * float(dist[(y, state)][key]) / 100
               for y in years)

print('\n=== Arm comparison (pooled) ===')
t1 = []
def add(label, op, nonop):
    print(f'{label:<28} operative {op:>8}   non-operative {nonop:>8}')
    t1.append({'characteristic': label, 'operative': op, 'non_operative': nonop})

add('Episodes, n', Ns, Nc)
add('Male, %', round(100*wsum('operative', 'male_pct')/Ns, 1),
    round(100*(wsum('cohort', 'male_pct') - wsum('operative', 'male_pct'))/Nc, 1))
# Three age groups; the full nine-band distribution is in age_specific.csv and Figure 1.
AGE_GROUPS = (('Age <65 years, %', AGE_KEYS[:6]),
              ('Age 65-79 years, %', ['age6574', 'age7579']),
              ('Age >=80 years, %', ['age80p']))
for lab, bands in AGE_GROUPS:
    Nso = sum(pool[b][0] for b in bands)
    Nco = sum(pool[b][2] for b in bands)
    add(lab, round(100*Nso/Ns, 1), round(100*Nco/Nc, 1))
losc = sum(int(res[y]['N_total'])*float(dist[(y, 'cohort')]['vwd_mean']) for y in years)
loss = sum(int(res[y]['N_operative'])*float(dist[(y, 'operative')]['vwd_mean']) for y in years)
add('Mean length of stay, days', round(loss/Ns, 1), round((losc-loss)/Nc, 1))

# Acute conditions associated with operative treatment, then chronic comorbidity.
if panel_in is not None:
    o_ = sum(int(panel_in[(y, 'Q2')]['G06.1']) for y in years)
    a_ = sum(int(panel_in[(y, 'Q1')]['G06.1']) for y in years)
    add('Intraspinal abscess, %', round(100*o_/Ns, 1), round(100*(a_-o_)/Nc, 1))
if exact_in is not None:
    for code, lab in (('G82', 'Paraparesis or paraplegia, %'), ('A41', 'Sepsis, %')):
        o_ = sum(exact_in[(y, 'operative', code, 'all')] for y in years)
        a_ = sum(exact_in[(y, 'cohort', code, 'all')] for y in years)
        add(lab, round(100*o_/Ns, 1), round(100*(a_-o_)/Nc, 1))
if panel_in is not None:
    o_ = sum(int(panel_in[(y, 'Q2')]['R57.2']) for y in years)
    a_ = sum(int(panel_in[(y, 'Q1')]['R57.2']) for y in years)
    add('Septic shock, %', round(100*o_/Ns, 1), round(100*(a_-o_)/Nc, 1))
if exact_in is not None:
    for code, lab in (('I50', 'Heart failure, %'), ('N18', 'Chronic kidney disease, %'),
                      ('E11', 'Type 2 diabetes, %')):
        o_ = sum(exact_in[(y, 'operative', code, 'all')] for y in years)
        a_ = sum(exact_in[(y, 'cohort', code, 'all')] for y in years)
        add(lab, round(100*o_/Ns, 1), round(100*(a_-o_)/Nc, 1))

o = sum(wsum('operative', k) for k in ['pccl4', 'pccl5', 'pccl6'])
c = sum(wsum('cohort', k) for k in ['pccl4', 'pccl5', 'pccl6'])
add('PCCL 4-6, %', round(100*o/Ns, 1), round(100*(c-o)/Nc, 1))
for lab, kS, kT in (('Intensive care, %', 'ICU_operative', 'ICU_total'),
                    ('Ventilation >96 h, %', 'vent96_operative', 'vent96_total')):
    S_, T_ = g(kS), g(kT)
    add(lab, round(100*S_/Ns, 1), round(100*(T_-S_)/Nc, 1))
add('In-hospital deaths, n', Ds, Dc)
add('In-hospital mortality, %', round(100*Ds/Ns, 2), round(100*Dc/Nc, 2))
write('arm_comparison.csv', t1)

# ------------------------------------------------------- 2. association measures
print('\n=== Crude and age-adjusted association by data year ===')
rows = []
year_est = []
for y in years:
    cells = [(s[2], s[3], s[4], s[5]) for s in strata if s[0] == y]
    m, lo, hi, se = mh_rr(cells)
    r = res[y]
    nso, dso = int(r['N_operative']), int(r['deaths_operative'])
    nco, dco = int(r['N_total'])-nso, int(r['deaths_total'])-dso
    c, clo, chi = katz_rr(dso, nso, dco, nco)
    year_est.append((math.log(m), se))
    print(f'{y}: crude RR {c:.3f} ({clo:.3f}-{chi:.3f})   '
          f'age-adjusted RR {m:.3f} ({lo:.3f}-{hi:.3f})')
    rows.append({'year': y, 'crude_rr': round(c, 6), 'crude_lo': round(clo, 6),
                 'crude_hi': round(chi, 6), 'mh_rr': round(m, 6), 'mh_lo': round(lo, 6),
                 'mh_hi': round(hi, 6)})

cells_all = [(s[2], s[3], s[4], s[5]) for s in strata]
m, lo, hi, _ = mh_rr(cells_all)
c, clo, chi = katz_rr(Ds, Ns, Dc, Nc)
print(f'pooled: crude RR {c:.3f} ({clo:.3f}-{chi:.3f})   '
      f'age-adjusted RR {m:.3f} ({lo:.3f}-{hi:.3f})  [stratified by year x age band]')
rows.append({'year': 'pooled', 'crude_rr': round(c, 6), 'crude_lo': round(clo, 6),
             'crude_hi': round(chi, 6), 'mh_rr': round(m, 6), 'mh_lo': round(lo, 6),
             'mh_hi': round(hi, 6)})
write('risk_ratios.csv', rows)

# heterogeneity across data years
wts = [1/se**2 for _, se in year_est]
theta = sum(w*t for w, (t, _) in zip(wts, year_est)) / sum(wts)
Q = sum(w*(t-theta)**2 for w, (t, _) in zip(wts, year_est))
df = len(year_est) - 1
I2 = max(0.0, 100*(Q-df)/Q)
pQ = chi2_sf(Q, df)
print(f'\nbetween-year heterogeneity: Q = {Q:.2f} on {df} df, P = {pQ:.2f}, I2 = {I2:.0f}%')

# ------------------------------------------------- 3. direct age standardisation
print('\n=== Direct age standardisation to the pooled cohort age distribution ===')
tot_by_band = {b: pool[b][0] + pool[b][2] for b in AGE_KEYS}
W = sum(tot_by_band.values())
std = {}
for arm, iN, iD in (('operative', 0, 1), ('non_operative', 2, 3)):
    rate = var = 0.0
    for b in AGE_KEYS:
        w = tot_by_band[b] / W
        n, d = pool[b][iN], pool[b][iD]
        if n <= 0:
            continue
        p = d / n
        rate += w * p
        var += w*w * p*(1-p) / n
    std[arm] = (100*rate, 100*math.sqrt(var))
diff = std['operative'][0] - std['non_operative'][0]
se_d = math.sqrt(std['operative'][1]**2 + std['non_operative'][1]**2)
print(f"operative     {std['operative'][0]:.2f}% (95% CI "
      f"{std['operative'][0]-Z*std['operative'][1]:.2f}-{std['operative'][0]+Z*std['operative'][1]:.2f})")
print(f"non-operative {std['non_operative'][0]:.2f}% (95% CI "
      f"{std['non_operative'][0]-Z*std['non_operative'][1]:.2f}-{std['non_operative'][0]+Z*std['non_operative'][1]:.2f})")
print(f'difference    {diff:.2f} percentage points (95% CI {diff-Z*se_d:.2f}-{diff+Z*se_d:.2f})')
write('standardised_mortality.csv', [
    {'quantity': 'directly standardised mortality, operative, %',
     'estimate': round(std['operative'][0], 2),
     'lo': round(std['operative'][0]-Z*std['operative'][1], 2),
     'hi': round(std['operative'][0]+Z*std['operative'][1], 2)},
    {'quantity': 'directly standardised mortality, non-operative, %',
     'estimate': round(std['non_operative'][0], 2),
     'lo': round(std['non_operative'][0]-Z*std['non_operative'][1], 2),
     'hi': round(std['non_operative'][0]+Z*std['non_operative'][1], 2)},
    {'quantity': 'standardised difference, percentage points',
     'estimate': round(diff, 2), 'lo': round(diff-Z*se_d, 2), 'hi': round(diff+Z*se_d, 2)},
    {'quantity': 'between-year heterogeneity of the age-adjusted RR, Cochran Q',
     'estimate': round(Q, 2), 'lo': '', 'hi': ''},
    {'quantity': 'between-year heterogeneity, I2, %', 'estimate': round(I2, 0),
     'lo': '', 'hi': ''},
    {'quantity': 'between-year heterogeneity, P value', 'estimate': round(pQ, 3),
     'lo': '', 'hi': ''},
])

# --------------------------------------------- 4. rounding-envelope sensitivity
print('\n=== Rounding-envelope sensitivity (published percentages are 2 dp) ===')
rng = random.Random(20260813)
draws = int(os.environ.get('SPONDYLO_DRAWS', '2000'))
vals = []
for _ in range(draws):
    s_, _c = build_strata(perturb=lambda p: p + rng.uniform(-0.005, 0.005))
    vals.append(mh_rr([(s[2], s[3], s[4], s[5]) for s in s_])[0])
vals.sort()
q = lambda f: vals[min(len(vals)-1, int(f*len(vals)))]
print(f'{draws} draws: median {vals[len(vals)//2]:.4f}, '
      f'2.5th-97.5th percentile {q(0.025):.4f}-{q(0.975):.4f}, '
      f'range {vals[0]:.4f}-{vals[-1]:.4f}')
write('rounding_envelope.csv', [{'draws': draws, 'point_estimate': round(m, 4),
                                 'median': round(vals[len(vals)//2], 4),
                                 'pct2_5': round(q(0.025), 4), 'pct97_5': round(q(0.975), 4),
                                 'min': round(vals[0], 4), 'max': round(vals[-1], 4)}])
if clipped:
    print('derived cells clipped at zero (rounding noise):', clipped)

# ------------------------------- 4b. distribution-residual sensitivity
print('\n=== Distribution-residual sensitivity ===')
base_dist = {k: dict(v) for k, v in dist.items()}


def mh_with(distmod):
    """Recompute the pooled Mantel-Haenszel estimate under a modified distribution set."""
    global dist
    dist = distmod
    s_, _c = build_strata()
    out = mh_rr([(x[2], x[3], x[4], x[5]) for x in s_])[0]
    dist = base_dist
    return out


bad = []
for (y_, state_), v in base_dist.items():
    tot_pct = sum(float(v[k]) for k in AGE_KEYS)
    if abs(tot_pct - 100) > 0.05:
        bad.append((y_, state_, tot_pct))
resid_rows = []
if bad:
    print('transcribed age distributions not summing to 100 (percentage points):')
    per_row = {}
    for y_, state_, tot_pct in bad:
        vals = []
        resid = 100 - tot_pct
        for b in AGE_KEYS:
            dm = {k: dict(v) for k, v in base_dist.items()}
            newv = float(dm[(y_, state_)][b]) + resid
            if newv < 0:
                continue
            dm[(y_, state_)][b] = str(newv)
            vals.append((mh_with(dm), b))
        vals.sort()
        per_row[(y_, state_)] = vals
        print(f'  {y_} {state_}: sum {tot_pct:.2f}; attributing the residual to any single '
              f'band gives RR {vals[0][0]:.4f} to {vals[-1][0]:.4f}')
        resid_rows.append({'year': y_, 'state': state_, 'published_sum_pct': round(tot_pct, 2),
                           'rr_min': round(vals[0][0], 4), 'rr_min_band': vals[0][1],
                           'rr_max': round(vals[-1][0], 4), 'rr_max_band': vals[-1][1]})
    joint = {}
    for label, idx in (('min', 0), ('max', -1)):
        dm = {k: dict(v) for k, v in base_dist.items()}
        for y_, state_, tot_pct in bad:
            b = per_row[(y_, state_)][idx][1]
            dm[(y_, state_)][b] = str(float(dm[(y_, state_)][b]) + (100 - tot_pct))
        joint[label] = mh_with(dm)
    print(f'  joint extremes across all affected rows: {joint["min"]:.4f} to {joint["max"]:.4f}')
    resid_rows.append({'year': 'joint', 'state': 'all affected rows', 'published_sum_pct': '',
                       'rr_min': round(joint['min'], 4), 'rr_min_band': '',
                       'rr_max': round(joint['max'], 4), 'rr_max_band': ''})
    write('distribution_residual.csv', resid_rows)
else:
    print('all transcribed age distributions sum to 100 within 0.05 percentage points')

# ----------------------------------------------------- 5. exposure sensitivity
print('\n=== Exposure-definition sensitivity (pooled) ===')
sens = []
_dbp = os.path.join(DATA, 'exposure_debridement.csv')
debr = {int(r['year']): r for r in rd('exposure_debridement.csv')} if os.path.exists(_dbp) else None
defs = [('core (primary)', 'N_operative', 'deaths_operative', 'complement'),
        ('wide (core + lumbar access)', 'N_operative_wide', 'deaths_operative_wide', 'complement'),
        ('instrumented only', 'N_instrumented', 'deaths_instrumented', 'core non-operative')]
if debr is not None:
    defs.append(('debridement-directed (excludes 5-839 and 5-83b)',
                 'N_debridement', 'deaths_debridement', 'core non-operative'))
for name, kN, kD, comparator in defs:
    if kN.endswith('debridement'):
        Ne = sum(int(debr[y][kN]) for y in years)
        De = sum(int(debr[y][kD]) for y in years)
    else:
        Ne, De = g(kN), g(kD)
    if comparator == 'complement':
        Nr, Dr = Nt - Ne, Dt - De
    else:
        Nr, Dr = Nc, Dc
    r_, l_, h_ = katz_rr(De, Ne, Dr, Nr)
    print(f'{name:<30} {De}/{Ne} = {100*De/Ne:.2f}%  vs {Dr}/{Nr} = {100*Dr/Nr:.2f}%  '
          f'crude RR {r_:.3f} ({l_:.3f}-{h_:.3f})')
    sens.append({'definition': name, 'comparator': comparator, 'N_exposed': Ne,
                 'deaths_exposed': De, 'mortality_pct': round(100*De/Ne, 2),
                 'N_comparator': Nr, 'deaths_comparator': Dr, 'rr': round(r_, 4),
                 'lo': round(l_, 4), 'hi': round(h_, 4)})
write('exposure_sensitivity.csv', sens)

# ------------------------------------- 6. indirect standardisation and E-values
print('\n=== Indirect standardisation and E-values (supplementary) ===')
smr_rows = []
obs_all = exp_all = 0.0
for y in years + ['pooled']:
    if y == 'pooled':
        obs, expd = obs_all, exp_all
    else:
        obs = expd = 0.0
        for _, band, Nso, Dso, Nco, Dco in [s for s in strata if s[0] == y]:
            obs += Dso
            expd += (Dco/Nco if Nco > 0 else 0.0) * Nso
        obs_all += obs; exp_all += expd
    o = obs
    smr = o/expd
    lo_b = o*(1 - 1/(9*o) - Z/(3*math.sqrt(o)))**3 / expd
    hi_b = (o+1)*(1 - 1/(9*(o+1)) + Z/(3*math.sqrt(o+1)))**3 / expd
    print(f'{y}: observed {obs:.0f}, expected {expd:.1f}, SMR {smr:.3f} ({lo_b:.3f}-{hi_b:.3f})')
    smr_rows.append({'year': y, 'observed': round(obs, 1), 'expected': round(expd, 1),
                     'smr': round(smr, 6), 'lo': round(lo_b, 6), 'hi': round(hi_b, 6)})
write('smr.csv', smr_rows)
print(f'E-value for the pooled age-adjusted RR {m:.3f}: {e_value(m):.2f}; '
      f'for the interval limit nearer the null {lo:.3f}: {e_value(lo):.2f}')
write('e_values.csv', [
    {'quantity': 'pooled age-adjusted risk ratio', 'estimate': round(m, 6),
     'e_value': round(e_value(m), 6)},
    {'quantity': 'interval limit nearer the null', 'estimate': round(lo, 6),
     'e_value': round(e_value(lo), 6)}])

# ------------------------------------------------------- 7. hospital structure
print('\n=== Hospital structure (suppressed residual cells excluded) ===')
def bed_group(cls):
    if cls in ('bis 49 Betten', '50 bis 99 Betten', '100 bis 149 Betten',
               '150 bis 199 Betten', '200 bis 249 Betten', '250 bis 299 Betten'):
        return '<300'
    if cls in ('300 bis 399 Betten', '400 bis 499 Betten', '500 bis 599 Betten'):
        return '300-599'
    if cls in ('600 bis 799 Betten', '800 bis 999 Betten'):
        return '600-999'
    if cls == '1000 Betten und mehr':
        return '1000+'
    return 'Rest'

agg = defaultdict(lambda: defaultdict(float))
rest = defaultdict(float)
for r in kv:
    q_, n_ = r['query'], float(r['faelle'])
    if r['suppressed_rest'] == 'True':
        rest[q_] += n_
        continue
    for dim in (('owner', r['traegerschaft']), ('beds', bed_group(r['bettengroessenklasse']))):
        agg[dim][q_] += n_
prov = []
for (kind, name), qs in sorted(agg.items()):
    N1, N2, N3, N4 = qs.get('Q1', 0), qs.get('Q2', 0), qs.get('Q3', 0), qs.get('Q4', 0)
    share = 100*N2/N1
    ms = 100*N4/N2 if N2 else float('nan')
    mc = 100*(N3-N4)/(N1-N2) if N1 > N2 else float('nan')
    print(f'{kind:>6} {name:<18} operative share {share:5.1f}%   '
          f'mortality operative {ms:5.2f}%   non-operative {mc:5.2f}%')
    prov.append({'dimension': kind, 'category': name, 'N_total': round(N1),
                 'N_operative': round(N2), 'operative_share_pct': round(share, 1),
                 'deaths_total': round(N3), 'deaths_operative': round(N4),
                 'mort_operative_pct': round(ms, 2), 'mort_nonoperative_pct': round(mc, 2)})
print('suppressed residual episodes excluded per query:', {k: int(v) for k, v in rest.items()})
write('hospital_structure.csv', prov)

# ------------------------------------------------- 7b. procedure-family counts
fam_path = os.path.join(DATA, 'procedure_family_counts.csv')
if os.path.exists(fam_path):
    print('\n=== Episodes identified by each operative OPS family ===')
    fam = list(csv.DictReader(open(fam_path)))
    by_fam = defaultdict(lambda: {'episodes': 0, 'codes': set()})
    for r in fam:
        f = by_fam[r['family']]
        f['episodes'] += int(r['episodes'])
        f['codes'].add(int(r['codes']))
    fam_rows = []
    for name in sorted(by_fam):
        f = by_fam[name]
        cmin, cmax = min(f['codes']), max(f['codes'])
        print(f"{name:>7}: {f['episodes']:6d} episodes  ({100*f['episodes']/Ns:5.1f}% of the operative arm)  "
              f"terminal codes {cmin}-{cmax}")
        fam_rows.append({'family': name, 'terminal_codes_min': cmin, 'terminal_codes_max': cmax,
                         'episodes_2019_2024': f['episodes'],
                         'pct_of_operative_arm': round(100*f['episodes']/Ns, 1)})
    total_fam = sum(r['episodes_2019_2024'] for r in fam_rows)
    print(f'sum over families {total_fam} against {Ns} operative episodes '
          f'(ratio {total_fam/Ns:.2f}); families are not mutually exclusive, an episode may '
          f'carry codes from several')
    fam_rows.append({'family': 'sum over families', 'terminal_codes_min': '', 'terminal_codes_max': '',
                     'episodes_2019_2024': total_fam,
                     'pct_of_operative_arm': round(100*total_fam/Ns, 1)})
    write('procedure_families.csv', fam_rows)

# ------------------------------- 7c. cohort-definition sensitivity (any diagnosis)
sec_path = os.path.join(DATA, 'secondary_diagnosis_totals.csv')
if os.path.exists(sec_path):
    print('\n=== Cohort definition: principal versus secondary diagnosis ===')
    sec = list(csv.DictReader(open(sec_path)))
    gs = lambda k: sum(int(r[k]) for r in sec)
    sN, sO = gs('N_total'), gs('N_operative')
    sD, sDO = gs('deaths_total'), gs('deaths_operative')
    ov_path = os.path.join(DATA, 'cohort_overlap_totals.csv')
    if os.path.exists(ov_path):
        ov = list(csv.DictReader(open(ov_path)))
        go = lambda k: sum(int(r[k]) for r in ov)
        oN, oO, oD, oDO = go('N_total'), go('N_operative'), go('deaths_total'), go('deaths_operative')
        spec = (('secondary diagnosis only (overlap removed)',
                 (sN-oN, sO-oO, sD-oD, sDO-oDO)),
                ('overlap (both principal and secondary)', (oN, oO, oD, oDO)),
                ('any diagnosis (union)', (Nt+sN-oN, Ns+sO-oO, Dt+sD-oD, Ds+sDO-oDO)))
    else:
        spec = (('any diagnosis (sum of the two, overlap not measured)',
                 (Nt+sN, Ns+sO, Dt+sD, Ds+sDO)),)
    rows = []
    for label, (tN, tO, tD, tDO) in ((
            ('principal diagnosis (primary cohort)', (Nt, Ns, Dt, Ds)),
            ('secondary diagnosis (overlap included)', (sN, sO, sD, sDO))) + spec):
        r_, l_, h_ = katz_rr(tDO, tO, tD-tDO, tN-tO)
        print(f'{label:<38} N={tN:6d}  operative {100*tO/tN:4.1f}%  mortality {100*tD/tN:5.2f}%  '
              f'operative {100*tDO/tO:5.2f}% vs non-operative {100*(tD-tDO)/(tN-tO):5.2f}%  '
              f'crude RR {r_:.3f} ({l_:.3f}-{h_:.3f})')
        rows.append({'cohort': label, 'N_total': tN, 'N_operative': tO,
                     'operative_share_pct': round(100*tO/tN, 1), 'deaths_total': tD,
                     'mortality_pct': round(100*tD/tN, 2),
                     'mort_operative_pct': round(100*tDO/tO, 2),
                     'mort_nonoperative_pct': round(100*(tD-tDO)/(tN-tO), 2),
                     'rr': round(r_, 4), 'lo': round(l_, 4), 'hi': round(h_, 4)})
    if os.path.exists(ov_path):
        print(f'the union removes the {oN} episodes carrying one M46.2-M46.5 code as principal '
              f'and another as secondary ({100*oN/Nt:.1f}% of the principal cohort)')
        per_year = [(int(r['year']),
                     int(r['N_total']) - int(q['N_total']),
                     int(r['N_operative']) - int(q['N_operative']),
                     int(r['deaths_total']) - int(q['deaths_total']),
                     int(r['deaths_operative']) - int(q['deaths_operative']))
                    for r, q in zip(sec, ov)]
        assert all(min(v[1:]) >= 0 for v in per_year), \
            'the overlap is not contained in the secondary-diagnosis cohort'
        print('secondary-only counts are non-negative in every data year')
    else:
        print('the any-diagnosis row treats the two cohorts as disjoint')

    # age-adjusted comparison within the secondary-diagnosis cohort
    sdist_path = os.path.join(DATA, 'secondary_diagnosis_distributions.csv')
    if os.path.exists(sdist_path):
        sdist = {(r['year'], r['state']): r for r in
                 csv.DictReader(open(sdist_path))}

        def sage(year, state, total):
            d = sdist[(year, state)]
            raw = [float(d[k])/100.0*total for k in AGE_KEYS]
            t = sum(raw)
            return [x*total/t for x in raw] if t > 0 else raw

        scells = []
        for r in sec:
            y = r['year']
            tN, tO = int(r['N_total']), int(r['N_operative'])
            tD, tDO = int(r['deaths_total']), int(r['deaths_operative'])
            ct, cs = sage(y, 'cohort', tN), sage(y, 'operative', tO)
            dt, ds = sage(y, 'deaths_total', tD), sage(y, 'deaths_operative', tDO)
            for i in range(len(AGE_KEYS)):
                No, Do = cs[i], ds[i]
                Nn, Dn = max(ct[i]-cs[i], 0.0), max(dt[i]-ds[i], 0.0)
                scells.append((No, min(Do, No), Nn, min(Dn, Nn)))
        sm, slo, shi, _ = mh_rr(scells)
        print(f'secondary-diagnosis cohort, age-adjusted RR {sm:.3f} ({slo:.3f}-{shi:.3f})')
        rows.append({'cohort': 'secondary diagnosis (overlap included), age-adjusted',
                     'N_total': sN, 'N_operative': sO, 'operative_share_pct': round(100*sO/sN, 1),
                     'deaths_total': sD, 'mortality_pct': round(100*sD/sN, 2),
                     'mort_operative_pct': round(100*sDO/sO, 2),
                     'mort_nonoperative_pct': round(100*(sD-sDO)/(sN-sO), 2),
                     'rr': round(sm, 4), 'lo': round(slo, 4), 'hi': round(shi, 4)})
    write('cohort_definition.csv', rows)

# ------------------------------------------- 7d. diagnosis-specificity analysis
spec_path = os.path.join(DATA, 'diagnosis_specificity.csv')
if os.path.exists(spec_path):
    print('\n=== Specific versus unspecified diagnosis codes ===')
    sp = list(csv.DictReader(open(spec_path)))
    gp = lambda k: sum(int(r[k]) for r in sp)
    spN, unN = gp('N_specific'), gp('N_unspecific')
    spD, unD = gp('deaths_specific'), gp('deaths_unspecific')
    assert spN + unN == Nt and spD + unD == Dt, 'specificity split does not reconcile'
    r_, l_, h_ = katz_rr(spD, spN, unD, unN)
    print(f'specific (M46.2, M46.3):    {spN:6d} episodes ({100*spN/Nt:4.1f}%), '
          f'mortality {100*spD/spN:.2f}%')
    print(f'unspecified (M46.4, M46.5): {unN:6d} episodes ({100*unN/Nt:4.1f}%), '
          f'mortality {100*unD/unN:.2f}%')
    print(f'mortality RR specific versus unspecified {r_:.3f} ({l_:.3f}-{h_:.3f})')
    write('diagnosis_specificity.csv', [
        {'code_group': 'specific (M46.2, M46.3)', 'episodes': spN,
         'share_of_cohort_pct': round(100*spN/Nt, 1), 'deaths': spD,
         'mortality_pct': round(100*spD/spN, 2), 'rr_vs_unspecified': round(r_, 4),
         'lo': round(l_, 4), 'hi': round(h_, 4)},
        {'code_group': 'unspecified (M46.4, M46.5)', 'episodes': unN,
         'share_of_cohort_pct': round(100*unN/Nt, 1), 'deaths': unD,
         'mortality_pct': round(100*unD/unN, 2), 'rr_vs_unspecified': 1,
         'lo': '', 'hi': ''}])

# ------------------------------- 7e. indication, comorbidity and discharge destination
panel_path = os.path.join(DATA, 'secondary_diagnosis_panel.csv')
exact_path = os.path.join(DATA, 'comorbidity_exact.csv')
disch_path = os.path.join(DATA, 'discharge_destination.csv')

if panel_in is not None and exact_in is not None:
    print('\n=== Coded clinical conditions and comorbidity (secondary-diagnosis register) ===')
    panel, ex = panel_in, exact_in

    def arm_counts(y):
        r = res[y]
        return (int(r['N_operative']), int(r['N_total']) - int(r['N_operative']),
                int(r['deaths_operative']),
                int(r['deaths_total']) - int(r['deaths_operative']))

    # (a) single ICD codes read directly from the register; these are exact
    #     because the code has no subcodes, subject to the suppression floor of 5.
    single = [('G06.1', 'Intraspinal / epidural abscess (G06.1)'),
              ('R57.2', 'Septic shock (R57.2)')]
    comorb_rows = []
    reg = {}
    for code, label in single:
        n1 = sum(int(panel[(y, 'Q2')][code]) for y in years)
        n0 = sum(int(panel[(y, 'Q1')][code]) - int(panel[(y, 'Q2')][code]) for y in years)
        d1 = sum(int(panel[(y, 'Q4')][code]) for y in years)
        d0 = sum(int(panel[(y, 'Q3')][code]) - int(panel[(y, 'Q4')][code]) for y in years)
        reg[code] = (n1, n0, d1, d0)
        pr, plo, phi = katz_rr(n1, Ns, n0, Nc)
        print(f'{label:<44} operative {100*n1/Ns:5.2f}%  non-operative {100*n0/Nc:5.2f}%  '
              f'prevalence ratio {pr:.2f} ({plo:.2f}-{phi:.2f})')
        comorb_rows.append({'condition': label, 'source': 'register (single code)',
                            'N_operative': n1, 'pct_operative': round(100*n1/Ns, 2),
                            'N_nonoperative': n0, 'pct_nonoperative': round(100*n0/Nc, 2),
                            'prevalence_ratio': round(pr, 3), 'lo': round(plo, 3),
                            'hi': round(phi, 3)})

    # (b) code groups counted by filter queries; the register cannot be summed
    #     across subcodes because an episode may carry several of them.
    groups = [('A41', 'Sepsis (A41.-)'), ('G82', 'Paraparesis / paraplegia (G82.-)'),
              ('E11', 'Type 2 diabetes (E11.-)'), ('I50', 'Heart failure (I50.-)'),
              ('N18', 'Chronic kidney disease (N18.-)')]
    for gcode, label in groups:
        n1 = sum(ex[(y, 'operative', gcode, 'all')] for y in years)
        n0 = sum(ex[(y, 'cohort', gcode, 'all')] - ex[(y, 'operative', gcode, 'all')]
                 for y in years)
        pr, plo, phi = katz_rr(n1, Ns, n0, Nc)
        print(f'{label:<44} operative {100*n1/Ns:5.2f}%  non-operative {100*n0/Nc:5.2f}%  '
              f'prevalence ratio {pr:.2f} ({plo:.2f}-{phi:.2f})')
        comorb_rows.append({'condition': label, 'source': 'filter query (code group)',
                            'N_operative': n1, 'pct_operative': round(100*n1/Ns, 2),
                            'N_nonoperative': n0, 'pct_nonoperative': round(100*n0/Nc, 2),
                            'prevalence_ratio': round(pr, 3), 'lo': round(plo, 3),
                            'hi': round(phi, 3)})

    # (c) coded secondary diagnosis entries per episode.
    #     The register's Nennungen column counts entries, not distinct codes: a code
    #     recorded in more than one department segment of the same episode is counted more
    #     than once. The Faelle column counts episodes per code, so summing it gives a
    #     distinct-code count, but only a lower bound, because the suppressed residual row
    #     reports episodes carrying at least one suppressed code rather than code-episode
    #     pairs. Both are reported; the arm contrast holds on either.
    nn1 = sum(int(panel[(y, 'Q2')]['nenn']) for y in years)
    nn0 = sum(int(panel[(y, 'Q1')]['nenn']) - int(panel[(y, 'Q2')]['nenn']) for y in years)
    print(f'{"Coded secondary diagnosis entries per episode":<44} operative {nn1/Ns:5.2f}   '
          f'non-operative {nn0/Nc:5.2f}')
    comorb_rows.append({'condition': 'Coded secondary diagnosis entries per episode (mean)',
                        'source': 'register (all codes)', 'N_operative': nn1,
                        'pct_operative': round(nn1/Ns, 2), 'N_nonoperative': nn0,
                        'pct_nonoperative': round(nn0/Nc, 2), 'prevalence_ratio': '',
                        'lo': '', 'hi': ''})
    if 'faelle' in panel[(years[0], 'Q1')]:
        f1 = sum(int(panel[(y, 'Q2')]['faelle']) for y in years)
        f0 = sum(int(panel[(y, 'Q1')]['faelle']) - int(panel[(y, 'Q2')]['faelle']) for y in years)
        print(f'{"Distinct codes per episode (lower bound)":<44} operative {f1/Ns:5.2f}   '
              f'non-operative {f0/Nc:5.2f}   entry inflation {nn1/f1:.4f} vs {nn0/f0:.4f}')
        comorb_rows.append({'condition': 'Distinct secondary diagnosis codes per episode '
                                         '(mean, lower bound)',
                            'source': 'register (all codes)', 'N_operative': f1,
                            'pct_operative': round(f1/Ns, 2), 'N_nonoperative': f0,
                            'pct_nonoperative': round(f0/Nc, 2), 'prevalence_ratio': '',
                            'lo': '', 'hi': ''})
    write('indication_comorbidity.csv', comorb_rows)

    # over-count check: register subcode sums against the filter-query counts
    oc_rows = []
    for gcode, label in groups:
        s = 0
        for y in years:
            s += int(panel[(y, 'Q1')][gcode])
        e = sum(ex[(y, 'cohort', gcode, 'all')] for y in years)
        oc_rows.append({'group': gcode, 'register_subcode_sum': s, 'filter_query_count': e,
                        'excess_pct': round(100*(s-e)/e, 1)})
    print('register subcode sums exceed the filter-query counts by',
          ', '.join(f"{r['group']} {r['excess_pct']}%" for r in oc_rows))
    write('register_subcode_overcount.csv', oc_rows)

    # (d) adjustment sets: each covariate is marginal to the others, because the
    #     interface returns marginal tables and not the joint distribution.
    print('\n=== Adjustment sets (each marginal to the others) ===')
    adj_rows = []

    def add(label, cells, note):
        r_, l_, h_, _ = mh_rr(cells)
        print(f'{label:<46} RR {r_:.3f} ({l_:.3f}-{h_:.3f})')
        adj_rows.append({'adjustment': label, 'rr': round(r_, 4), 'lo': round(l_, 4),
                         'hi': round(h_, 4), 'note': note})
        return r_

    crude = katz_rr(Ds, Ns, Dc, Nc)
    print(f'{"none (crude)":<46} RR {crude[0]:.3f} ({crude[1]:.3f}-{crude[2]:.3f})')
    adj_rows.append({'adjustment': 'none (crude)', 'rr': round(crude[0], 4),
                     'lo': round(crude[1], 4), 'hi': round(crude[2], 4),
                     'note': 'Katz logarithmic interval'})
    yr_cells = [(int(res[y]['N_operative']), int(res[y]['deaths_operative']),
                 int(res[y]['N_total']) - int(res[y]['N_operative']),
                 int(res[y]['deaths_total']) - int(res[y]['deaths_operative']))
                for y in years]
    add('data year', yr_cells, 'six strata')
    add('data year and age band', [(s[2], s[3], s[4], s[5]) for s in strata],
        'primary analysis, 54 strata')

    def two_level(getter, name):
        cells = []
        for y in years:
            No_, Nn_, Do_, Dn_ = arm_counts(y)
            n1, n0, d1, d0 = getter(y)
            cells.append((n1, d1, n0, d0))
            cells.append((No_ - n1, Do_ - d1, Nn_ - n0, Dn_ - d0))
        return add(name, cells, '12 strata (condition x year)')

    def absc(y):
        return (int(panel[(y, 'Q2')]['G06.1']),
                int(panel[(y, 'Q1')]['G06.1']) - int(panel[(y, 'Q2')]['G06.1']),
                int(panel[(y, 'Q4')]['G06.1']),
                int(panel[(y, 'Q3')]['G06.1']) - int(panel[(y, 'Q4')]['G06.1']))

    def seps(y):
        return (ex[(y, 'operative', 'A41', 'all')],
                ex[(y, 'cohort', 'A41', 'all')] - ex[(y, 'operative', 'A41', 'all')],
                ex[(y, 'operative', 'A41', 'death')],
                ex[(y, 'cohort', 'A41', 'death')] - ex[(y, 'operative', 'A41', 'death')])

    two_level(absc, 'data year and intraspinal abscess (G06.1)')
    two_level(seps, 'data year and sepsis (A41.-)')
    write('adjustment_sets.csv', adj_rows)

    # (e) how much confounding each measured indication can generate
    print('\n=== Bounding factors of the measured indications (VanderWeele-Ding) ===')
    bf_rows = []
    for label, (n1, n0, d1, d0) in (
            ('intraspinal abscess (G06.1)', reg['G06.1']),
            ('septic shock (R57.2)', reg['R57.2']),
            ('sepsis (A41.-)', (sum(ex[(y, 'operative', 'A41', 'all')] for y in years),
                                sum(ex[(y, 'cohort', 'A41', 'all')]
                                    - ex[(y, 'operative', 'A41', 'all')] for y in years),
                                sum(ex[(y, 'operative', 'A41', 'death')] for y in years),
                                sum(ex[(y, 'cohort', 'A41', 'death')]
                                    - ex[(y, 'operative', 'A41', 'death')] for y in years)))):
        rr_eu = (n1/Ns) / (n0/Nc)
        rr_ud = ((d1 + d0)/(n1 + n0)) / ((Dt - d1 - d0)/(Nt - n1 - n0))
        bf = rr_eu*rr_ud/(rr_eu + rr_ud - 1)
        print(f'{label:<32} exposure association {rr_eu:.2f}, mortality association '
              f'{rr_ud:.2f}, bounding factor {bf:.2f}')
        bf_rows.append({'condition': label, 'rr_exposure_confounder': round(rr_eu, 3),
                        'rr_confounder_outcome': round(rr_ud, 3), 'bounding_factor': round(bf, 3)})
    write('bounding_factors.csv', bf_rows)

if os.path.exists(disch_path):
    print('\n=== Discharge destination ===')
    dd = {}
    for r in csv.DictReader(open(disch_path)):
        dd[(int(r['year']), r['arm'], r['destination'])] = r['n']
    dest = [('06+08', 'Transfer to another acute hospital'),
            ('09', 'Discharge to rehabilitation'),
            ('10', 'Discharge to a nursing facility'),
            ('11', 'Discharge to a hospice')]
    drows = []
    for code, label in dest:
        vals = [dd[(y, a, code)] for y in years for a in ('cohort', 'operative')]
        sup = any(v == 'SUP' for v in vals)
        n1 = sum(int(dd[(y, 'operative', code)]) for y in years
                 if dd[(y, 'operative', code)] != 'SUP')
        nc = sum(int(dd[(y, 'cohort', code)]) for y in years
                 if dd[(y, 'cohort', code)] != 'SUP')
        n0 = nc - n1
        if sup:
            print(f'{label:<38} cohort {100*nc/Nt:5.2f}%  (arm-specific counts suppressed '
                  f'in at least one year)')
            drows.append({'destination': label, 'N_operative': '', 'pct_operative': '',
                          'N_nonoperative': '', 'pct_nonoperative': '', 'N_cohort': nc,
                          'pct_cohort': round(100*nc/Nt, 2), 'ratio': ''})
            continue
        print(f'{label:<38} operative {100*n1/Ns:5.2f}%  non-operative {100*n0/Nc:5.2f}%  '
              f'ratio {(n1/Ns)/(n0/Nc):.2f}')
        drows.append({'destination': label, 'N_operative': n1,
                      'pct_operative': round(100*n1/Ns, 2), 'N_nonoperative': n0,
                      'pct_nonoperative': round(100*n0/Nc, 2), 'N_cohort': nc,
                      'pct_cohort': round(100*nc/Nt, 2),
                      'ratio': round((n1/Ns)/(n0/Nc), 3)})
    write('discharge_destination.csv', drows)

    # transfer sensitivity: in-hospital death is censored by transfer out, and
    # transfer is more common in the non-operative arm.
    def rr_after_transfer(p1, p0):
        cells = []
        for y in years:
            No_ = int(res[y]['N_operative']); Nn_ = int(res[y]['N_total']) - No_
            Do_ = int(res[y]['deaths_operative'])
            Dn_ = int(res[y]['deaths_total']) - Do_
            t1 = int(dd[(y, 'operative', '06+08')])
            t0 = int(dd[(y, 'cohort', '06+08')]) - t1
            cells.append((No_, Do_ + p1*t1, Nn_, Dn_ + p0*t0))
        return mh_rr(cells)[0]

    print('\n=== Transfer-out sensitivity (year-stratified risk ratio) ===')
    trows = []
    for p in (0.0, 0.02, 0.05, 0.10, 0.20):
        r_ = rr_after_transfer(p, p)
        print(f'post-transfer mortality {100*p:5.1f}% in both arms -> RR {r_:.3f}')
        trows.append({'scenario': f'{100*p:.0f}% in both arms', 'rr': round(r_, 4)})
    lo_p, hi_p = 0.0, 1.0
    for _ in range(60):
        mid = (lo_p + hi_p)/2
        if rr_after_transfer(mid, mid) > 1.0:
            lo_p = mid
        else:
            hi_p = mid
    print(f'equal post-transfer mortality that would null the association: {100*lo_p:.1f}%')
    trows.append({'scenario': 'equal rate that nulls the association',
                  'rr': round(100*lo_p, 1)})
    lo_p2, hi_p2 = 0.0, 1.0
    for _ in range(60):
        mid = (lo_p2 + hi_p2)/2
        if rr_after_transfer(0.0, mid) > 1.0:
            lo_p2 = mid
        else:
            hi_p2 = mid
    print(f'post-transfer mortality in the non-operative arm alone that would null it: '
          f'{100*lo_p2:.1f}%')
    trows.append({'scenario': 'non-operative arm only, rate that nulls the association',
                  'rr': round(100*lo_p2, 1)})
    write('transfer_sensitivity.csv', trows)

# ------------------------------------------------------------ 8. annual summary
annual = []
for y in years:
    r = res[y]
    N, O = int(r['N_total']), int(r['N_operative'])
    D, DO = int(r['deaths_total']), int(r['deaths_operative'])
    d_ = dist[(y, 'cohort')]
    annual.append({'year': y, 'episodes': N, 'operative': O,
                   'operative_share_pct': round(100*O/N, 1),
                   'male_pct': d_['male_pct'], 'age80plus_pct': d_['age80p'],
                   'mean_los_days': d_['vwd_mean'],
                   'pccl4_6_pct': round(float(d_['pccl4'])+float(d_['pccl5'])+float(d_['pccl6']), 1),
                   'icu_pct': round(100*int(r['ICU_total'])/N, 1),
                   'vent96_pct': round(100*int(r['vent96_total'])/N, 1),
                   'deaths': D, 'mortality_pct': round(100*D/N, 2),
                   'deaths_operative': DO, 'mort_operative_pct': round(100*DO/O, 2),
                   'mort_nonoperative_pct': round(100*(D-DO)/(N-O), 2)})
write('annual_summary.csv', annual)
print('\ndone')
