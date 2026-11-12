#!/usr/bin/env python3
"""Figures 1 and 2 and eFigure 1 of the article, from the derived results.

Reads results/figure1_data.csv, results/figure2_data.csv and results/risk_ratios.csv and
writes vector PDF and 300 dpi PNG to figures/. Requires matplotlib. Deterministic.
"""
import csv, os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator

HERE = Path(__file__).resolve().parent
RES = Path(os.environ.get('SPONDYLO_RESULTS', HERE.parent / 'results'))
OUT = Path(os.environ.get('SPONDYLO_FIGURES', HERE.parent / 'figures'))
OUT.mkdir(exist_ok=True)

# Nord palette: nord10 for operative, nord3 for nonoperative, nord0 for text.
BLUE, GREY, INK, RULE = '#5E81AC', '#4C566A', '#2E3440', '#D8DEE9'
plt.rcParams.update({
    'font.family': ['Liberation Sans', 'DejaVu Sans'], 'font.size': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5,
    'ytick.major.size': 2.5, 'axes.edgecolor': INK, 'axes.labelcolor': INK,
    'xtick.color': INK, 'ytick.color': INK, 'text.color': INK,
    'pdf.fonttype': 42, 'savefig.dpi': 300,
})
CM = 1 / 2.54
WIDTH = 8.5 * CM  # single column


def rows(name):
    return list(csv.DictReader(open(RES / name, encoding='utf-8')))


def tidy(ax):
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    ax.spines['bottom'].set_position(('outward', 3))


def save(fig, stem):
    fig.savefig(OUT / f'{stem}.pdf', bbox_inches='tight', pad_inches=0.02,
                metadata={'CreationDate': None})
    fig.savefig(OUT / f'{stem}.png', bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


def errs(vals, lo, hi):
    return [[v - l for v, l in zip(vals, lo)], [h - v for v, h in zip(vals, hi)]]


# ------------------------------------------------------------------ Figure 1
f1 = [r for r in rows('figure1_data.csv') if r['age_group'] != '<30']
labels = [r['age_group'].replace('80+', '≥80').replace('-', '–') for r in f1]
x = list(range(len(f1)))
share = [float(r['operative_share_pct']) for r in f1]
m1 = [float(r['op_mort_pct']) for r in f1]; l1 = [float(r['op_lo']) for r in f1]; h1 = [float(r['op_hi']) for r in f1]
m0 = [float(r['nonop_mort_pct']) for r in f1]; l0 = [float(r['nonop_lo']) for r in f1]; h0 = [float(r['nonop_hi']) for r in f1]

fig, (a, b) = plt.subplots(2, 1, figsize=(WIDTH, 9.2 * CM), sharex=True,
                           gridspec_kw={'height_ratios': [1, 1.45], 'hspace': 0.18})
a.plot(x, share, color=BLUE, lw=0.8, marker='o', ms=3.2, mfc='white', mec=BLUE, mew=0.8)
a.set_ylim(30, 60); a.yaxis.set_major_locator(FixedLocator([30, 40, 50, 60]))
a.set_ylabel('Admissions with\nsurgery, %')
a.text(-0.21, 1.0, 'A', transform=a.transAxes, fontsize=7.5, va='top')
tidy(a); a.tick_params(axis='x', length=0)
a.spines['bottom'].set_visible(False)

off = 0.13
b.errorbar([i - off for i in x], m1, yerr=errs(m1, l1, h1), fmt='o', ms=3.0, mfc=BLUE, mec=BLUE,
           ecolor=BLUE, elinewidth=0.7, capsize=0)
b.errorbar([i + off for i in x], m0, yerr=errs(m0, l0, h0), fmt='s', ms=2.8, mfc='white', mec=GREY,
           mew=0.8, ecolor=GREY, elinewidth=0.7, capsize=0)
b.set_ylim(0, 16); b.yaxis.set_major_locator(FixedLocator([0, 4, 8, 12, 16]))
b.set_ylabel('In-hospital mortality, %')
b.set_xticks(x); b.set_xticklabels(labels)
b.set_xlabel('Age, years')
b.text(-0.21, 1.0, 'B', transform=b.transAxes, fontsize=7.5, va='top')
# direct labels at the oldest age group
b.text(x[-1] - off, h1[-1] + 0.45, 'With surgery', color=BLUE, fontsize=6.5, ha='center', va='bottom')
b.text(x[-1] + off, l0[-1] - 0.45, 'Without\nsurgery', color=GREY, fontsize=6.5, ha='center', va='top')
tidy(b)
save(fig, 'article_figure1_age')

# ------------------------------------------------------------------ Figure 2
f2 = rows('figure2_data.csv')
blab = ['<300', '300–599', '600–999', '≥1000']
xb = list(range(len(f2)))
s = [float(r['operative_share_pct']) for r in f2]
sl = [float(r['share_lo']) for r in f2]; sh = [float(r['share_hi']) for r in f2]
bm1 = [float(r['op_mort_pct']) for r in f2]; bl1 = [float(r['op_lo']) for r in f2]; bh1 = [float(r['op_hi']) for r in f2]
bm0 = [float(r['nonop_mort_pct']) for r in f2]; bl0 = [float(r['nonop_lo']) for r in f2]; bh0 = [float(r['nonop_hi']) for r in f2]

fig, (a, b) = plt.subplots(2, 1, figsize=(WIDTH, 8.6 * CM), sharex=True,
                           gridspec_kw={'height_ratios': [1, 1.2], 'hspace': 0.2})
a.errorbar(xb, s, yerr=errs(s, sl, sh), fmt='o', ms=3.2, mfc='white', mec=BLUE, mew=0.8,
           ecolor=BLUE, elinewidth=0.7, capsize=0)
a.plot(xb, s, color=BLUE, lw=0.6)
a.set_ylim(20, 70); a.yaxis.set_major_locator(FixedLocator([20, 30, 40, 50, 60, 70]))
a.set_ylabel('Admissions with\nsurgery, %')
a.text(-0.21, 1.0, 'A', transform=a.transAxes, fontsize=7.5, va='top')
tidy(a); a.tick_params(axis='x', length=0); a.spines['bottom'].set_visible(False)
a.set_xlim(-0.5, len(xb) - 0.5)

off = 0.11
b.errorbar([i - off for i in xb], bm1, yerr=errs(bm1, bl1, bh1), fmt='o', ms=3.0, mfc=BLUE, mec=BLUE,
           ecolor=BLUE, elinewidth=0.7, capsize=0)
b.errorbar([i + off for i in xb], bm0, yerr=errs(bm0, bl0, bh0), fmt='s', ms=2.8, mfc='white', mec=GREY,
           mew=0.8, ecolor=GREY, elinewidth=0.7, capsize=0)
b.set_ylim(0, 10); b.yaxis.set_major_locator(FixedLocator([0, 2, 4, 6, 8, 10]))
b.set_ylabel('In-hospital mortality, %')
b.set_xticks(xb); b.set_xticklabels(blab)
b.set_xlabel('Hospital size, beds')
b.text(-0.21, 1.0, 'B', transform=b.transAxes, fontsize=7.5, va='top')
b.text(xb[0] - off, bl1[0] - 0.3, 'With\nsurgery', color=BLUE, fontsize=6.5, ha='center', va='top')
b.text(xb[0] + off, bh0[0] + 0.3, 'Without\nsurgery', color=GREY, fontsize=6.5, ha='center', va='bottom')
tidy(b)
save(fig, 'article_figure2_hospital_size')

# ------------------------------------------------------------------ eFigure 1: annual risk ratios
rr = rows('risk_ratios.csv')
yrs = [r['year'] for r in rr]
fig, ax = plt.subplots(figsize=(WIDTH, 6.2 * CM))
yy = list(range(len(rr)))[::-1]
for y_, r in zip(yy, rr):
    pooled = r['year'] == 'pooled'
    ax.plot([float(r['mh_lo']), float(r['mh_hi'])], [y_, y_], color=BLUE, lw=0.8 if pooled else 0.7)
    ax.plot(float(r['mh_rr']), y_, marker='D' if pooled else 'o', ms=3.6 if pooled else 3.0,
            color=BLUE, mfc=BLUE)
    ax.plot(float(r['crude_rr']), y_ + 0.22, marker='s', ms=2.6, mfc='white', mec=GREY, mew=0.7, ls='none')
ax.axvline(1.0, color=RULE, lw=0.7, zorder=0)
ax.set_xscale('log'); ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_major_locator(FixedLocator([0.8, 1.0, 1.25, 1.5, 1.75]))
ax.set_xticklabels(['0.8', '1.0', '1.25', '1.5', '1.75'])
ax.set_xlim(0.78, 1.8)
ax.set_yticks(yy); ax.set_yticklabels([('Pooled' if y == 'pooled' else y) for y in yrs])
ax.set_xlabel('Risk ratio, surgery vs no surgery')
tidy(ax)
save(fig, 'article_efigure1_annual_risk_ratios')
