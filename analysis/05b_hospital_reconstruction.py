#!/usr/bin/env python3
"""Cell-level reconstruction of the hospital-structure tables (bed size x ownership).

The DatenBrowser hospital distribution ("Krankenhaus-Verteilung") hides some cells and
reports them only as a combined "Rest" row. Hiding is not restricted to small counts
(in 2020 the cohort query hid three large-hospital cells holding more than 240
admissions), and it differs between the four queries of a year (Q1 cohort, Q2
operative, Q3 deaths, Q4 operative deaths). Pooled tables built from the displayed
cells are therefore internally inconsistent: a cell can be present in the numerator
query and absent from the denominator query.

This script works at the level of year x bed-size class x ownership (36 cells per
year). For each query with a Rest row it uses the number of hidden rows and their
total, both transcribed in data/hospital_hidden_rows.csv, and identifies:
  - mandatory hidden cells: absent from the query but logically non-zero, because
    the cell is displayed with a positive count in a query that it must contain
    (for example present in Q2 but absent from Q1);
  - optional hidden cells: absent from the query and possibly non-zero;
  - rows outside the 36-cell grid (for example an unclassified category), which
    exist because some years hide rows although all 36 grid cells are displayed.
Each draw selects the hidden rows, distributes the Rest total over them (each
hidden row >= 1), and accepts the draw only if every cell satisfies
Q2 <= Q1, Q3 <= Q1, Q4 <= Q2, Q4 <= Q3 and Q3 - Q4 <= Q1 - Q2. Mass assigned to rows
outside the grid cannot be attributed to a hospital category and is excluded from
the hospital analyses. 2000 accepted draws with the remainder spread uniformly (medians and percentiles)
and 2000 extreme draws with the remainder concentrated (min and max), fixed seed. Standard library only.

Inputs: data/hospital_tables.csv, data/hospital_hidden_rows.csv (SPONDYLO_DATA).
Outputs (results/): hospital_reconstruction_draws.csv, hospital_reconstruction_summary.csv,
hospital_reconstruction_cells.csv (identification of hidden cells, for audit).
"""
import csv, math, os, random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get('SPONDYLO_RESULTS', os.path.join(HERE, '..', 'results'))
DATA = os.environ.get('SPONDYLO_DATA', os.path.join(HERE, '..', 'data'))
SEED, NDRAWS, MAXTRY = 20261112, 2000, 400000
Z = 1.959963984540054
Q = ('Q1', 'Q2', 'Q3', 'Q4')
YEARS = ['2019', '2020', '2021', '2022', '2023', '2024']


def bed_group(cls):
    if cls in ('bis 49 Betten', '50 bis 99 Betten', '100 bis 149 Betten', '150 bis 199 Betten',
               '200 bis 249 Betten', '250 bis 299 Betten'):
        return '<300'
    if cls in ('300 bis 399 Betten', '400 bis 499 Betten', '500 bis 599 Betten'):
        return '300-599'
    if cls in ('600 bis 799 Betten', '800 bis 999 Betten'):
        return '600-999'
    if cls == '1000 Betten und mehr':
        return '1000+'
    raise ValueError(cls)


BEDS = ['<300', '300-599', '600-999', '1000+']
OWN = ['Oeffentlich', 'Freigemeinnuetzig', 'Privat']


def mh_rr(cells):
    num = den = P = 0.0
    for Ns, Ds, Nc, Dc in cells:
        T = Ns + Nc
        if T <= 0 or Ns <= 0 or Nc <= 0:
            continue
        num += Ds * Nc / T
        den += Dc * Ns / T
        P += ((Ds + Dc) * Ns * Nc - Ds * Dc * T) / T**2
    rr = num / den
    se = math.sqrt(P / (num * den))
    return rr, rr * math.exp(-Z * se), rr * math.exp(Z * se)


# ------------------------------------------------------------- inputs
kv = list(csv.DictReader(open(os.path.join(DATA, 'hospital_tables.csv'), encoding='utf-8')))
GRID = sorted({(r['bettengroessenklasse'], r['traegerschaft']) for r in kv if r['suppressed_rest'] == 'False'})
assert len(GRID) == 36, len(GRID)
shown = {y: {c: {} for c in GRID} for y in YEARS}
rest_total = {y: {} for y in YEARS}
for r in kv:
    if r['suppressed_rest'] == 'True':
        rest_total[r['year']][r['query']] = int(r['faelle'])
    else:
        shown[r['year']][(r['bettengroessenklasse'], r['traegerschaft'])][r['query']] = int(r['faelle'])

rest_rows = {y: {} for y in YEARS}
for r in csv.DictReader(open(os.path.join(DATA, 'hospital_hidden_rows.csv'), encoding='utf-8')):
    y, q = r['year'], r['query']
    k, tot = int(r['hidden_rows']), int(r['hidden_admissions'])
    assert rest_total[y].get(q) == tot, (y, q, tot, rest_total[y].get(q))
    rest_rows[y][q] = k
for y in YEARS:  # every Rest total has a row count
    assert set(rest_rows[y]) == set(rest_total[y]), y

# ------------------------------------------------------------- hidden-cell identification
PARENT = {'Q2': ['Q1'], 'Q3': ['Q1'], 'Q4': ['Q2', 'Q3']}   # q is contained in its parents
CHILD = {'Q1': ['Q2', 'Q3', 'Q4'], 'Q2': ['Q4'], 'Q3': ['Q4'], 'Q4': []}
ident = []
plan = {y: {} for y in YEARS}
for y in YEARS:
    for q in rest_total[y]:
        missing = [c for c in GRID if q not in shown[y][c]]
        mandatory = [c for c in missing if any(shown[y][c].get(ch, 0) > 0 for ch in CHILD[q])]
        optional = []
        for c in missing:
            if c in mandatory:
                continue
            # a cell can be non-zero in q only if no displayed parent is zero
            par = [shown[y][c].get(p) for p in PARENT.get(q, [])]
            if any(v == 0 for v in par):
                continue
            optional.append(c)
        k, R = rest_rows[y][q], rest_total[y][q]
        assert len(mandatory) <= k, (y, q, mandatory, k)
        plan[y][q] = (k, R, mandatory, optional)
        for c in missing:
            ident.append({'year': y, 'query': q, 'bed_class': c[0], 'ownership': c[1],
                          'status': 'mandatory hidden' if c in mandatory else ('optional' if c in optional else 'zero'),
                          'rest_rows': k, 'rest_total': R})

# ------------------------------------------------------------- sampling
rng = random.Random(SEED)


def coherent(v):
    n1, n2, d, d1 = v['Q1'], v['Q2'], v['Q3'], v['Q4']
    return n2 <= n1 and d <= n1 and d1 <= n2 and d1 <= d and (d - d1) <= (n1 - n2)


def draw_year(y, vertex=False):
    vals = {c: {q: shown[y][c].get(q, 0) for q in Q} for c in GRID}
    outside = {q: 0 for q in Q}
    for q, (k, R, mand, opt) in plan[y].items():
        free = k - len(mand)
        pool = opt + ['__outside__'] * free      # each free row may lie outside the grid
        pick = list(mand) + rng.sample(pool, free) if free else list(mand)
        lbs, ubs = [], []
        for c in pick:
            if c == '__outside__':
                lbs.append(1); ubs.append(R)
                continue
            sv = shown[y][c]
            n1, n2, dd, d1 = (sv.get(x) for x in Q)
            if q == 'Q1':
                lb = max(1, (n2 or 0) + max((dd or 0) - (d1 or 0), 0), dd or 0)
                ub = R
            elif q == 'Q2':
                lb = max(1, d1 or 0)
                ub = (n1 - max((dd or 0) - (d1 or 0), 0)) if n1 is not None else R
            elif q == 'Q3':
                lb = max(1, d1 or 0)
                ub = n1 if n1 is not None else R
            else:
                lb = 1
                ub = min(x for x in (n2, dd, R) if x is not None)
            lbs.append(lb); ubs.append(ub)
        rem = R - sum(lbs)
        if rem < 0 or any(l > u for l, u in zip(lbs, ubs)):
            return None
        alloc = list(lbs)
        open_ = [i for i in range(k) if alloc[i] < ubs[i]]
        if vertex:   # extreme allocation: as much of the remainder as possible on one row
            order = open_[:]
            rng.shuffle(order)
            for i in order:
                add_ = min(rem, ubs[i] - alloc[i])
                alloc[i] += add_
                rem -= add_
            if rem:
                return None
            open_ = []
        for _ in range(rem):
            if not open_:
                return None
            i = rng.choice(open_)
            alloc[i] += 1
            if alloc[i] >= ubs[i]:
                open_.remove(i)
        for c, a in zip(pick, alloc):
            if c == '__outside__':
                outside[q] += a
            else:
                vals[c][q] += a
    for c in GRID:
        if not coherent(vals[c]):
            return None
    return vals, outside


# Constraints apply within a year, so years are sampled independently and
# combined draw by draw.
per_year, acc_rate = {}, {}
for y in YEARS:
    got, tries = [], 0
    while len(got) < NDRAWS:
        tries += 1
        if tries > MAXTRY:
            raise RuntimeError(f'{y}: only {len(got)} accepted draws after {MAXTRY} attempts')
        r_ = draw_year(y)
        if r_ is not None:
            got.append(r_)
    per_year[y], acc_rate[y] = got, len(got) / tries
draws = [{y: per_year[y][i] for y in YEARS} for i in range(NDRAWS)]

# Extreme (vertex) reconstructions: each Rest total is concentrated on as few hidden
# rows as the bounds allow. Their range bounds the estimates more widely than the
# uniform draws above and is reported as the feasible range.
vert = {}
for y in YEARS:
    got, tries = [], 0
    while len(got) < NDRAWS:
        tries += 1
        if tries > MAXTRY:
            raise RuntimeError(f'{y}: vertex sampling failed')
        r_ = draw_year(y, vertex=True)
        if r_ is not None:
            got.append(r_)
    vert[y] = got
vdraws = [{y: vert[y][i] for y in YEARS} for i in range(NDRAWS)]


def summarise(d):
    out = {}
    tot = defaultdict(lambda: [0, 0, 0, 0])
    yb, yo, bb, oo = [], [], [], []
    by_yb = defaultdict(lambda: [0, 0, 0, 0]); by_yo = defaultdict(lambda: [0, 0, 0, 0])
    for y in YEARS:
        vals, outside = d[y]
        for c in GRID:
            v = vals[c]
            for key in (('beds', bed_group(c[0])), ('owner', c[1])):
                t = tot[key]
                for i, q in enumerate(Q):
                    t[i] += v[q]
            for store, key in ((by_yb, (y, bed_group(c[0]))), (by_yo, (y, c[1]))):
                for i, q in enumerate(Q):
                    store[key][i] += v[q]
        for q in Q:
            out[f'outside_{q}'] = out.get(f'outside_{q}', 0) + outside[q]
    def cells(src, keys):
        return [(src[k][1], src[k][3], src[k][0] - src[k][1], src[k][2] - src[k][3]) for k in keys]
    out['rr_year_beds'] = mh_rr(cells(by_yb, [(y, b) for y in YEARS for b in BEDS]))
    out['rr_beds'] = mh_rr(cells(tot, [('beds', b) for b in BEDS]))
    out['rr_year_owner'] = mh_rr(cells(by_yo, [(y, o) for y in YEARS for o in OWN]))
    out['rr_owner'] = mh_rr(cells(tot, [('owner', o) for o in OWN]))
    a = sum(tot[('beds', b)][3] for b in BEDS); n1 = sum(tot[('beds', b)][1] for b in BEDS)
    b_ = sum(tot[('beds', b)][2] - tot[('beds', b)][3] for b in BEDS)
    n0 = sum(tot[('beds', b)][0] - tot[('beds', b)][1] for b in BEDS)
    out['rr_crude_grid'] = (a / n1) / (b_ / n0)
    for key, t in tot.items():
        out[('N1',) + key], out[('N2',) + key], out[('D',) + key], out[('D1',) + key] = t
    return out


S = [summarise(d) for d in draws]
SV = [summarise(d) for d in vdraws]


def pct(vals, p):
    v = sorted(vals)
    i = (len(v) - 1) * p
    lo, hi = math.floor(i), math.ceil(i)
    return v[lo] + (v[hi] - v[lo]) * (i - lo)


rows = []
def add(name, vals, vvals=None):
    allv = list(vals) + list(vvals or [])
    rows.append({'quantity': name, 'median': round(pct(vals, .5), 4), 'p2_5': round(pct(vals, .025), 4),
                 'p97_5': round(pct(vals, .975), 4), 'min': round(min(allv), 4), 'max': round(max(allv), 4)})

for key in ('rr_year_beds', 'rr_beds', 'rr_year_owner', 'rr_owner'):
    add(key, [s[key][0] for s in S], [s[key][0] for s in SV])
    # 95% CI of the draw whose estimate is the median
    med = pct([s[key][0] for s in S], .5)
    j = min(range(len(S)), key=lambda i: abs(S[i][key][0] - med))
    rows.append({'quantity': key + '_ci_at_median_draw', 'median': round(S[j][key][0], 4),
                 'p2_5': round(S[j][key][1], 4), 'p97_5': round(S[j][key][2], 4), 'min': '', 'max': ''})
add('rr_crude_grid', [s['rr_crude_grid'] for s in S])
for dim, cats in (('beds', BEDS), ('owner', OWN)):
    for c in cats:
        k = (dim, c)
        add(f'{dim}|{c}|admissions', [s[('N1',) + k] for s in S])
        add(f'{dim}|{c}|operative', [s[('N2',) + k] for s in S])
        add(f'{dim}|{c}|deaths', [s[('D',) + k] for s in S])
        add(f'{dim}|{c}|deaths_operative', [s[('D1',) + k] for s in S])
        add(f'{dim}|{c}|operative_share_pct', [100 * s[('N2',) + k] / s[('N1',) + k] for s in S],
            [100 * s[('N2',) + k] / s[('N1',) + k] for s in SV])
        add(f'{dim}|{c}|mort_operative_pct', [100 * s[('D1',) + k] / s[('N2',) + k] for s in S])
        add(f'{dim}|{c}|mort_nonoperative_pct',
            [100 * (s[('D',) + k] - s[('D1',) + k]) / (s[('N1',) + k] - s[('N2',) + k]) for s in S])
for q in Q:
    add(f'outside_grid_{q}', [s[f'outside_{q}'] for s in S])
gridN = [sum(s[('N1', 'beds', b)] for b in BEDS) for s in S]
add('grid_admissions', gridN)
add('share_admissions_lt300_pct', [100 * s[('N1', 'beds', '<300')] / g for s, g in zip(S, gridN)])
add('share_operations_ge600_pct', [100 * (s[('N2', 'beds', '600-999')] + s[('N2', 'beds', '1000+')]) /
                                   sum(s[('N2', 'beds', b)] for b in BEDS) for s in S])
rows.append({'quantity': 'accepted_draws', 'median': len(draws), 'p2_5': '', 'p97_5': '', 'min': '', 'max': ''})
for y in YEARS:
    rows.append({'quantity': f'acceptance_rate_{y}', 'median': round(acc_rate[y], 4), 'p2_5': '', 'p97_5': '', 'min': '', 'max': ''})


def write(name, rr):
    with open(os.path.join(RES, name), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rr[0].keys()))
        w.writeheader(); w.writerows(rr)


write('hospital_reconstruction_summary.csv', rows)
write('hospital_reconstruction_cells.csv', ident)
write('hospital_reconstruction_draws.csv',
      [{'draw': i, 'rr_year_beds': round(s['rr_year_beds'][0], 6), 'rr_beds': round(s['rr_beds'][0], 6),
        'rr_year_owner': round(s['rr_year_owner'][0], 6), 'rr_owner': round(s['rr_owner'][0], 6)}
       for i, s in enumerate(S)])
for r in rows:
    print(r)
