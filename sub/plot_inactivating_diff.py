#!/usr/bin/env python3
"""
Histograms of (my inactivating insertions - published inactivating insertions)
for all genes in the screen and for hits only (published p < 1e-5).
One panel per screen (PrP, CD59).
"""

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PUBLISHED = "/pl/active/ShenLab_PL/haploid_genetics/1-s2.0-S2211124715005525-mmc2.xlsx"
CD59_RES  = ("/scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens"
             "/CD59_vs_JScontrol_comparative_screen/CD59_vs_JScontrol_results.tsv")
PRP_RES   = ("/scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens"
             "/PrP_vs_JScontrol_comparative_screen/PrP_vs_JScontrol_results.tsv")
OUT_FIG   = "inactivating_diff.pdf"

HIT_PVAL  = 1e-5
BIN_WIDTH = 2   # insertions per bin

# ── Load and merge ────────────────────────────────────────────────────────────
pub = pd.read_excel(PUBLISHED, sheet_name='Sheet2', header=None, skiprows=3)
pub.columns = ['gene','prp_intron','prp_inact','prp_pval',
               'cd59_intron','cd59_inact','cd59_pval',
               'ctrl_intron','ctrl_inact']
pub = pub.dropna(subset=['gene']).copy()
pub['gene']      = pub['gene'].astype(str).str.strip()
pub['prp_pval']  = pd.to_numeric(pub['prp_pval'],  errors='coerce')
pub['cd59_pval'] = pd.to_numeric(pub['cd59_pval'], errors='coerce')
pub['prp_inact'] = pd.to_numeric(pub['prp_inact'], errors='coerce')
pub['cd59_inact']= pd.to_numeric(pub['cd59_inact'],errors='coerce')

my_prp  = pd.read_csv(PRP_RES,  sep='\t')
my_cd59 = pd.read_csv(CD59_RES, sep='\t')

screens = [
    ('PrP',  my_prp,  'prp_inact',  'prp_pval'),
    ('CD59', my_cd59, 'cd59_inact', 'cd59_pval'),
]

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, (label, my_df, inact_col, pval_col) in zip(axes, screens):

    merged = my_df.merge(
        pub[['gene', inact_col, pval_col]], on='gene', how='inner'
    )
    merged['diff'] = merged['inactivating_exp'] - merged[inact_col]
    hits = merged[merged[pval_col] < HIT_PVAL]

    lo = int(np.floor(merged['diff'].min() / BIN_WIDTH) * BIN_WIDTH)
    hi = int(np.ceil( merged['diff'].max() / BIN_WIDTH) * BIN_WIDTH)
    bins = np.arange(lo - BIN_WIDTH, hi + 2 * BIN_WIDTH, BIN_WIDTH)

    # All genes — density so it's comparable to the tiny hit set
    ax.hist(merged['diff'], bins=bins, density=True,
            color='steelblue', alpha=0.55, label=f'All genes (n={len(merged):,})',
            edgecolor='steelblue', linewidth=0.3)

    # Hits
    ax.hist(hits['diff'], bins=bins, density=True,
            color='tomato', alpha=0.75, label=f'Published hits (n={len(hits)}, p<1e-5)',
            edgecolor='darkred', linewidth=0.5)

    ax.axvline(0, color='black', linewidth=1, linestyle='--', alpha=0.6)

    ax.set_xlabel('My inactivating insertions − published inactivating insertions',
                  fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title(f'{label} Screen', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25, axis='y')

    # Median annotation for each group
    for diff_series, color, va in [
        (merged['diff'], 'steelblue', 'top'),
        (hits['diff'],   'tomato',    'bottom'),
    ]:
        med = diff_series.median()
        ax.axvline(med, color=color, linewidth=1.2, linestyle=':', alpha=0.9)
        ax.text(med, ax.get_ylim()[1] * 0.95, f'  median={med:.0f}',
                color=color, fontsize=8, va='top')

fig.suptitle(
    'Distribution of inactivating insertion count difference\n'
    '(this study − Davis et al. 2015)',
    fontsize=13, fontweight='bold', y=1.02,
)
plt.tight_layout()
plt.savefig(OUT_FIG, dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_FIG}")

# ── Summary table ─────────────────────────────────────────────────────────────
for label, my_df, inact_col, pval_col in screens:
    merged = my_df.merge(pub[['gene', inact_col, pval_col]], on='gene', how='inner')
    merged['diff'] = merged['inactivating_exp'] - merged[inact_col]
    hits = merged[merged[pval_col] < HIT_PVAL]
    print(f"\n{label} — all genes (n={len(merged):,}): "
          f"median={merged['diff'].median():.1f}, "
          f"mean={merged['diff'].mean():.2f}, "
          f"range=[{merged['diff'].min():.0f}, {merged['diff'].max():.0f}]")
    print(f"{label} — hits (n={len(hits)}): "
          f"median={hits['diff'].median():.1f}, "
          f"mean={hits['diff'].mean():.2f}, "
          f"range=[{hits['diff'].min():.0f}, {hits['diff'].max():.0f}]")
    print(hits[['gene','inactivating_exp', inact_col, 'diff']].sort_values('diff').to_string())
