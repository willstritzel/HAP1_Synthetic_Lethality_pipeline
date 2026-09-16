#!/usr/bin/env python3
"""
Cumulative rank overlap: pipeline results vs Davis et al. 2015 published results.

For each rank threshold N (1..30):
  y = |my top-N genes ∩ published top-N genes| / N * 100

Published genes are ranked by screen-specific p-value (PrP column for PrP plot,
CD59 column for CD59 plot).
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PUBLISHED    = "/pl/active/ShenLab_PL/haploid_genetics/1-s2.0-S2211124715005525-mmc2.xlsx"
CD59_RESULTS = ("/scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens"
                "/CD59_vs_JScontrol_comparative_screen/CD59_vs_JScontrol_results.tsv")
PRP_RESULTS  = ("/scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens"
                "/PrP_vs_JScontrol_comparative_screen/PrP_vs_JScontrol_results.tsv")
OUT_FIG      = "rank_overlap.pdf"

MAX_RANK = 30

# ── Load published results ────────────────────────────────────────────────────
pub = pd.read_excel(PUBLISHED, sheet_name='Sheet2', header=None, skiprows=3)
pub.columns = [
    'gene',
    'prp_intron', 'prp_inact',  'prp_pval',
    'cd59_intron', 'cd59_inact', 'cd59_pval',
    'ctrl_intron', 'ctrl_inact',
]
pub = pub.dropna(subset=['gene']).copy()
pub['gene']     = pub['gene'].astype(str).str.strip()
pub['prp_pval'] = pd.to_numeric(pub['prp_pval'],  errors='coerce')
pub['cd59_pval']= pd.to_numeric(pub['cd59_pval'], errors='coerce')

# Rank separately for each screen
pub_prp_ranked  = pub.sort_values('prp_pval' )['gene'].tolist()
pub_cd59_ranked = pub.sort_values('cd59_pval')['gene'].tolist()

# ── Load my results ───────────────────────────────────────────────────────────
my_prp  = pd.read_csv(PRP_RESULTS,  sep='\t').sort_values('p_value')
my_cd59 = pd.read_csv(CD59_RESULTS, sep='\t').sort_values('p_value')

my_prp_ranked  = my_prp['gene'].tolist()
my_cd59_ranked = my_cd59['gene'].tolist()

# ── Compute cumulative overlap ────────────────────────────────────────────────
def cumulative_overlap(my_ranked, pub_ranked, max_rank):
    out = []
    for n in range(1, max_rank + 1):
        my_set  = set(my_ranked[:n])
        pub_set = set(pub_ranked[:n])
        out.append(len(my_set & pub_set) / n * 100)
    return out

ranks     = list(range(1, MAX_RANK + 1))
prp_pcts  = cumulative_overlap(my_prp_ranked,  pub_prp_ranked,  MAX_RANK)
cd59_pcts = cumulative_overlap(my_cd59_ranked, pub_cd59_ranked, MAX_RANK)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

for ax, pcts, title in zip(
    axes,
    [prp_pcts, cd59_pcts],
    ['PrP Screen', 'CD59 Screen'],
):
    ax.plot(ranks, pcts, marker='o', linewidth=2, markersize=5, color='steelblue')
    ax.set_xlabel('Rank threshold (top N genes)', fontsize=12)
    ax.set_ylabel('% overlap with published top-N', fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlim(0.5, MAX_RANK + 0.5)
    ax.set_ylim(0, 105)
    ax.set_xticks(range(5, MAX_RANK + 1, 5))
    ax.set_xticks(range(1, MAX_RANK + 1), minor=True)
    ax.grid(True, alpha=0.3)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(decimals=0))

    # Annotate each point with its value
    for n, pct in zip(ranks, pcts):
        if pct > 0:
            ax.annotate(
                f'{pct:.0f}%',
                xy=(n, pct),
                xytext=(0, 6),
                textcoords='offset points',
                ha='center',
                fontsize=6.5,
                color='steelblue',
            )

fig.suptitle(
    'Cumulative rank overlap vs Davis et al. 2015',
    fontsize=14, fontweight='bold', y=1.02,
)
plt.tight_layout()
plt.savefig(OUT_FIG, dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_FIG}")

# ── Print top-30 gene lists for inspection ────────────────────────────────────
print("\nPrP — my top 30:", my_prp_ranked[:30])
print("PrP — published top 30:", pub_prp_ranked[:30])
print("\nCD59 — my top 30:", my_cd59_ranked[:30])
print("CD59 — published top 30:", pub_cd59_ranked[:30])
