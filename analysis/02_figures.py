#!/usr/bin/env python3
"""Figures for the treatment-selection analysis.

Style: Tufte-minimal. No titles inside the figure, thin axes, greyscale with a
single muted accent, open markers, direct labelling, vector PDF plus 300 dpi PNG.
"""
import csv, math, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator, NullFormatter

RES = os.environ.get('SPONDYLO_RESULTS', os.path.join(os.path.dirname(__file__), '..', 'results'))
FIG = os.environ.get('SPONDYLO_FIGURES', os.path.join(os.path.dirname(__file__), '..', 'figures'))
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 8, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
    'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
})
GREY = '#4d4d4d'
ACCENT = '#8c6bb1'
rd = lambda n: list(csv.DictReader(open(os.path.join(RES, n))))
save = lambda fig, stem: (fig.savefig(os.path.join(FIG, stem + '.pdf')),
                          fig.savefig(os.path.join(FIG, stem + '.png'), dpi=300),
                          plt.close(fig))


def wilson(k, n):
    if n == 0:
        return float('nan'), float('nan')
    p, z = k/n, 1.96
    d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return c-h, c+h


# ---- Figure 1: operative share by age (A) and mortality by age and arm (B) ----
age = [r for r in rd('age_specific.csv') if r['band'] != 'age_u30']
lbl = [r['label'].replace('-', '–').replace('80+', '≥80') for r in age]
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.4, 4.6), sharex=True,
                               gridspec_kw={'height_ratios': [1, 1.25]})
x = range(len(age))
share = [float(r['operative_share_pct']) for r in age]
ax1.plot(x, share, marker='o', ms=3.6, mfc='white', mec=ACCENT, mew=1.0,
         color=ACCENT, lw=0.8)
ax1.set_ylabel('Operative share, %')
ax1.set_ylim(25, 60)
ax1.set_yticks([30, 40, 50, 60])

for i, r in enumerate(age):
    for (N, D, col, dx) in ((int(r['N_operative']), int(r['D_operative']), ACCENT, -0.12),
                            (int(r['N_nonoperative']), int(r['D_nonoperative']), GREY, 0.12)):
        lo, hi = wilson(D, N)
        ax2.plot([i+dx, i+dx], [100*lo, 100*hi], color=col, lw=0.7)
        ax2.plot(i+dx, 100*D/N, marker='o', ms=3.4, mfc='white', mec=col, mew=0.9)
ax2.set_xticks(list(x)); ax2.set_xticklabels(lbl, fontsize=7.5)
ax2.set_xlabel('Age, years'); ax2.set_ylabel('In-hospital mortality, %')
ax2.set_ylim(0, 16)
ax2.text(0.03, 0.93, 'operative', color=ACCENT, transform=ax2.transAxes, fontsize=8)
ax2.text(0.03, 0.84, 'non-operative', color=GREY, transform=ax2.transAxes, fontsize=8)
for ax, tag in ((ax1, 'A'), (ax2, 'B')):
    ax.text(-0.13, 1.02, tag, transform=ax.transAxes, fontsize=8.5, va='top')
fig.tight_layout(h_pad=1.2)
save(fig, 'figure1_selection_and_mortality_by_age')

# ------------------- Figure 2: operative share by hospital size ---------------
prov = [r for r in rd('hospital_structure.csv') if r['dimension'] == 'beds']
order = ['<300', '300-599', '600-999', '1000+']
prov.sort(key=lambda r: order.index(r['category']))
fig, ax = plt.subplots(figsize=(3.4, 2.4))
xs = range(len(prov))
share = [float(r['operative_share_pct']) for r in prov]
ax.plot(xs, share, marker='o', ms=4, mfc='white', mec=ACCENT, mew=1.0, color=ACCENT, lw=0.8)
for i, s_ in zip(xs, share):
    ax.annotate(f'{s_:.0f}%', (i, s_), textcoords='offset points', xytext=(0, 6),
                ha='center', fontsize=7.5, color=ACCENT)
ax.set_xticks(list(xs))
ax.set_xticklabels(['<300', '300–599', '600–999', '≥1000'])
ax.set_xlabel('Hospital size, beds'); ax.set_ylabel('Operative share, %')
ax.set_ylim(0, 75)
fig.tight_layout()
save(fig, 'figureS1_operative_share_by_hospital_size')

# ---------------- Figure 2: crude and age-adjusted risk ratios by year --------
rr = rd('risk_ratios.csv')
fig, ax = plt.subplots(figsize=(3.4, 2.6))
ys = []
for i, r in enumerate(rr):
    y = len(rr) - 1 - i
    ys.append((y, r['year']))
    pooled = r['year'] == 'pooled'
    if not pooled:
        ax.plot(float(r['crude_rr']), y+0.16, marker='o', ms=3.2, mfc='white',
                mec=GREY, mew=0.9)
    m, lo, hi = float(r['mh_rr']), float(r['mh_lo']), float(r['mh_hi'])
    yy = y if pooled else y-0.16
    ax.plot([lo, hi], [yy, yy], color=ACCENT, lw=0.8)
    ax.plot(m, yy, marker='D' if pooled else 'o', ms=4 if pooled else 3.2,
            mfc=ACCENT if pooled else 'white', mec=ACCENT, mew=0.9)
ax.axvline(1.0, color='#bbbbbb', lw=0.6, zorder=0)
ax.set_yticks([y for y, _ in ys]); ax.set_yticklabels([n for _, n in ys])
ax.set_xscale('log')
ax.xaxis.set_major_locator(FixedLocator([0.8, 1.0, 1.25, 1.6]))
ax.xaxis.set_major_formatter(FixedFormatter(['0.8', '1.0', '1.25', '1.6']))
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlim(0.75, 1.8)
ax.set_xlabel('Risk ratio, operative vs non-operative (log scale)')
ax.text(0.02, 0.12, 'age-adjusted', color=ACCENT, transform=ax.transAxes, fontsize=8)
ax.text(0.02, 0.04, 'crude', color=GREY, transform=ax.transAxes, fontsize=8)
fig.tight_layout()
save(fig, 'figure2_risk_ratios_by_year')
# --------- Figure 3: clinical selection profile (prevalence ratios by arm) ----
ind_path = os.path.join(RES, 'indication_comorbidity.csv')
adj_path = os.path.join(RES, 'adjustment_sets.csv')
SHORT = {'Intraspinal / epidural abscess (G06.1)': 'Intraspinal abscess',
         'Septic shock (R57.2)': 'Septic shock',
         'Sepsis (A41.-)': 'Sepsis',
         'Paraparesis / paraplegia (G82.-)': 'Paraparesis or paraplegia',
         'Type 2 diabetes (E11.-)': 'Type 2 diabetes',
         'Heart failure (I50.-)': 'Heart failure',
         'Chronic kidney disease (N18.-)': 'Chronic kidney disease'}
if os.path.exists(ind_path):
    ind = [r for r in rd('indication_comorbidity.csv') if r['prevalence_ratio'] != '']
    ind = sorted(ind, key=lambda r: float(r['prevalence_ratio']))
    fig, ax = plt.subplots(figsize=(4.4, 2.9))
    for i, r in enumerate(ind):
        m, lo, hi = (float(r['prevalence_ratio']), float(r['lo']), float(r['hi']))
        ax.plot([lo, hi], [i, i], color=ACCENT, lw=0.8)
        ax.plot(m, i, marker='o', ms=4, mfc=ACCENT if m > 1 else 'white', mec=ACCENT, mew=1.0)
    ax.axvline(1.0, color='#bbbbbb', lw=0.6, zorder=0)
    ax.set_yticks(range(len(ind)))
    ax.set_yticklabels([SHORT[r['condition']] for r in ind])
    ax.set_ylim(-0.7, len(ind) - 0.3)
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(FixedLocator([0.7, 1.0, 1.5, 2.0, 3.0, 4.0]))
    ax.xaxis.set_major_formatter(FixedFormatter(['0.7', '1.0', '1.5', '2.0', '3.0', '4.0']))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(0.62, 4.6)
    ax.set_xlabel('Prevalence ratio, operative vs non-operative (log scale)')
    fig.tight_layout()
    save(fig, 'figure3_clinical_selection_profile')

# ---- Supplementary Figure S2: risk ratio under alternative stratification ----
if os.path.exists(adj_path):
    ORDER = ['none (crude)', 'data year', 'data year and age band',
             'data year and intraspinal abscess (G06.1)', 'data year and sepsis (A41.-)']
    NAME = {'none (crude)': 'None', 'data year': 'Year',
            'data year and age band': 'Year + age',
            'data year and intraspinal abscess (G06.1)': 'Year + abscess',
            'data year and sepsis (A41.-)': 'Year + sepsis'}
    adj = {r['adjustment']: r for r in rd('adjustment_sets.csv')}
    fig, ax = plt.subplots(figsize=(3.9, 2.4))
    for i, k in enumerate(ORDER):
        y = len(ORDER) - 1 - i
        r = adj[k]
        m, lo, hi = float(r['rr']), float(r['lo']), float(r['hi'])
        ax.plot([lo, hi], [y, y], color=ACCENT, lw=0.8)
        ax.plot(m, y, marker='o', ms=3.6, mfc='white', mec=ACCENT, mew=1.0)
    ax.axvline(1.0, color='#bbbbbb', lw=0.6, zorder=0)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels([NAME[k] for k in reversed(ORDER)])
    ax.set_ylim(-0.6, len(ORDER) - 0.4)
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(FixedLocator([0.85, 1.0, 1.2, 1.45]))
    ax.xaxis.set_major_formatter(FixedFormatter(['0.85', '1.0', '1.2', '1.45']))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(0.82, 1.5)
    ax.set_xlabel('Risk ratio, operative vs non-operative (log scale)')
    fig.tight_layout()
    save(fig, 'figureS2_stratification_sets')

print('figures written to', FIG)
