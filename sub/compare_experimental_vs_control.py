#!/usr/bin/env python3
"""
Compare inactivating insertion frequencies between experimental and control
haploid genetic screens, following Davis et al. 2015 (Cell Reports).

Inactivating insertions = intronic sense + all exonic (sense or antisense).
Fisher's exact test (one-sided: greater) tests whether each gene's share of
genome-wide inactivating insertions is significantly higher in the experimental
condition than in the unselected control.

Usage:
    python3 compare_experimental_vs_control.py \
        --intron-exp  <sense_vs_antisense_counts.csv> \
        --intron-ctrl <sense_vs_antisense_counts.csv> \
        --exon-exp    <exon_counts.tsv> \
        --exon-ctrl   <exon_counts.tsv> \
        --out         <results.tsv> \
        [--pthresh 1e-5]
"""
import argparse
import numpy as np
import pandas as pd
import scipy.stats
from statsmodels.stats.multitest import multipletests


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--intron-exp',  required=True,
                   help='Intronic sense_vs_antisense_counts.csv for experimental')
    p.add_argument('--intron-ctrl', required=True,
                   help='Intronic sense_vs_antisense_counts.csv for control')
    p.add_argument('--exon-exp',    required=True,
                   help='Exonic counts TSV for experimental (from count_exon_insertions.py)')
    p.add_argument('--exon-ctrl',   required=True,
                   help='Exonic counts TSV for control')
    p.add_argument('--out',         required=True, help='Output TSV path')
    p.add_argument('--pthresh',     type=float, default=1e-5,
                   help='P-value threshold to flag hits (default: 1e-5)')
    return p.parse_args()


def load_intron(path):
    df = pd.read_csv(path, sep='\t', usecols=['gene', 'sense', 'antisense'])
    df = df.rename(columns={'sense': 'intron_sense', 'antisense': 'intron_antisense'})
    return df


def load_exon(path):
    df = pd.read_csv(path, sep='\t', usecols=['gene', 'exon_insertions'])
    return df


def fisher_greater(i_exp, total_exp, i_ctrl, total_ctrl):
    """One-sided Fisher exact test: is i_exp / total_exp > i_ctrl / total_ctrl?"""
    table = [[i_exp,  total_exp  - i_exp],
             [i_ctrl, total_ctrl - i_ctrl]]
    _, p = scipy.stats.fisher_exact(table, alternative='greater')
    return p


def main():
    args = parse_args()

    # Load intronic counts
    exp_intron  = load_intron(args.intron_exp)
    ctrl_intron = load_intron(args.intron_ctrl)

    # Load exonic counts (fill missing genes with 0)
    exp_exon  = load_exon(args.exon_exp)
    ctrl_exon = load_exon(args.exon_ctrl)

    # Build per-condition frame: outer join intron + exon counts
    exp = exp_intron.merge(exp_exon, on='gene', how='outer').fillna(0)
    ctrl = ctrl_intron.merge(ctrl_exon, on='gene', how='outer').fillna(0)

    for col in ['intron_sense', 'intron_antisense', 'exon_insertions']:
        exp[col]  = exp[col].astype(int)
        ctrl[col] = ctrl[col].astype(int)

    # Merge experimental and control on gene (outer join)
    merged = exp.merge(ctrl, on='gene', how='outer', suffixes=('_exp', '_ctrl')).fillna(0)
    for col in merged.columns[1:]:
        merged[col] = merged[col].astype(int)

    # Inactivating insertions = intronic sense + all exonic
    merged['inactivating_exp']  = merged['intron_sense_exp']  + merged['exon_insertions_exp']
    merged['inactivating_ctrl'] = merged['intron_sense_ctrl'] + merged['exon_insertions_ctrl']

    # Genome-wide totals
    total_exp  = int(merged['inactivating_exp'].sum())
    total_ctrl = int(merged['inactivating_ctrl'].sum())
    print(f"Total inactivating insertions — experimental: {total_exp:,}, control: {total_ctrl:,}",
          flush=True)

    # Per-gene Fisher exact test (one-sided: enrichment in experimental)
    merged['p_value'] = merged.apply(
        lambda r: fisher_greater(
            int(r['inactivating_exp']),  total_exp,
            int(r['inactivating_ctrl']), total_ctrl
        ),
        axis=1
    )

    # Benjamini-Hochberg FDR across all genes
    _, fdr, _, _ = multipletests(merged['p_value'], method='fdr_bh')
    merged['fdr'] = fdr

    # Flag hits
    merged['hit'] = merged['p_value'] < args.pthresh

    # Sort by p_value ascending, then by inactivating_exp descending
    merged = merged.sort_values(['p_value', 'inactivating_exp'],
                                ascending=[True, False]).reset_index(drop=True)

    # Output columns
    out_cols = [
        'gene',
        'intron_sense_exp', 'intron_antisense_exp', 'exon_insertions_exp', 'inactivating_exp',
        'intron_sense_ctrl', 'intron_antisense_ctrl', 'exon_insertions_ctrl', 'inactivating_ctrl',
        'p_value', 'fdr', 'hit',
    ]
    merged[out_cols].to_csv(args.out, sep='\t', index=False)

    n_hits = merged['hit'].sum()
    print(f"Results written to {args.out}")
    print(f"Hits at p < {args.pthresh}: {n_hits} genes")


if __name__ == '__main__':
    main()
