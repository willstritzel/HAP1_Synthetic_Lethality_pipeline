#!/usr/bin/env python3
"""
plot_insertion_map.py — Visualize retroviral insertion site distributions
across gene bodies from HAP1 synthetic lethality screens.

Usage examples:
  # Single gene:
  python sub/plot_insertion_map.py --gene STX18

  # Top 5 genes from a p-value CSV:
  python sub/plot_insertion_map.py \
      --input /path/to/sense_vs_antisense_counts.csv --top 5

  # With custom BWT files or bin size:
  python sub/plot_insertion_map.py --gene STX18 \
      --bwt /path/to/replicate_*/replicate.mapped_unique.bwt \
      --bin-size 1000 --output-dir /scratch/alpine/wist9668/figures

Outputs (single-gene mode):
  {GENE}_insertion_map.pdf

Outputs (multi-gene mode, three files):
  {screen}_top{N}_individual.pdf          — one page per gene
  {screen}_top{N}_aligned_5prime.pdf      — stacked panels aligned by 5' end
  {screen}_top{N}_aligned_antisense_peak.pdf — stacked panels aligned by peak antisense bin
"""

import argparse
import bisect
import csv
import glob
import os
import subprocess
import sys

import matplotlib.ticker as ticker

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages

# ─── Chromosome name mappings ──────────────────────────────────────────────────
# NC (RefSeq) → chr name → CM (RefSeq complete genome) accession
NC_TO_CHR = {
    "NC_000001.11": "chr1",  "NC_000002.12": "chr2",  "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",  "NC_000005.10": "chr5",  "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",  "NC_000008.11": "chr8",  "NC_000009.12": "chr9",
    "NC_000010.11": "chr10", "NC_000011.10": "chr11", "NC_000012.12": "chr12",
    "NC_000013.11": "chr13", "NC_000014.9":  "chr14", "NC_000015.10": "chr15",
    "NC_000016.10": "chr16", "NC_000017.11": "chr17", "NC_000018.10": "chr18",
    "NC_000019.10": "chr19", "NC_000020.11": "chr20", "NC_000021.9":  "chr21",
    "NC_000022.11": "chr22", "NC_000023.11": "chrX",  "NC_000024.10": "chrY",
    "NC_012920.1":  "chrM",
}

CHR_TO_CM = {
    "chr1":  "CM000663.2", "chr2":  "CM000664.2", "chr3":  "CM000665.2",
    "chr4":  "CM000666.2", "chr5":  "CM000667.2", "chr6":  "CM000668.2",
    "chr7":  "CM000669.2", "chr8":  "CM000670.2", "chr9":  "CM000671.2",
    "chr10": "CM000672.2", "chr11": "CM000673.2", "chr12": "CM000674.2",
    "chr13": "CM000675.2", "chr14": "CM000676.2", "chr15": "CM000677.2",
    "chr16": "CM000678.2", "chr17": "CM000679.2", "chr18": "CM000680.2",
    "chr19": "CM000681.2", "chr20": "CM000682.2", "chr21": "CM000683.2",
    "chr22": "CM000684.2", "chrX":  "CM000685.2", "chrY":  "CM000686.2",
    "chrM":  "J01415.2",
}

NC_TO_CM = {nc: CHR_TO_CM[ch] for nc, ch in NC_TO_CHR.items() if ch in CHR_TO_CM}

# ─── Default file paths ────────────────────────────────────────────────────────
_DATA_BASE = (
    "/scratch/alpine/wist9668/data/analyzed_data/"
    "synthetic_lethal_screens/ControlData-HAP1_output"
)
DEFAULT_BWT_PATHS = [
    f"{_DATA_BASE}/replicate_{r}/replicate.mapped_unique.bwt"
    for r in range(1, 5)
]
DEFAULT_EXONS   = "/scratch/alpine/wist9668/reference/annotation/exons.bed"
DEFAULT_INTRONS = "/scratch/alpine/wist9668/reference/annotation/GRCh38_introns_CM_mapped.bed"
DEFAULT_GTF     = ("/scratch/alpine/wist9668/reference/annotation/"
                   "GCF_000001405.40_GRCh38.p14_genomic.gtf")

# ─── Colours ──────────────────────────────────────────────────────────────────
SENSE_COLOR    = "#2166ac"   # blue
ANTISENSE_COLOR = "#d7191c"  # red
EXON_FACE      = "#2c7bb6"
EXON_EDGE      = "#1a4f7a"


# ══════════════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════════════

def load_gene_structure(gene_name, exons_bed, introns_bed):
    """
    Return (chrom_cm, gene_start, gene_end, strand, exon_list, nc_chrom).
    Exon coords are sourced from exons.bed (NC accessions); the chromosome is
    mapped to CM accessions used by the BWT alignment files.
    exon_list: sorted list of (start, end) tuples, deduplicated.
    nc_chrom: the NC accession for the gene's chromosome (needed for GTF lookup).
    """
    exons = set()
    strand = None
    nc_chrom = None

    with open(exons_bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue
            if fields[3] != gene_name:
                continue
            exons.add((int(fields[1]), int(fields[2])))
            strand = fields[5]
            nc_chrom = fields[0]

    if not exons:
        # Fall back: try to find the gene in the CM introns BED
        with open(introns_bed) as fh:
            for line in fh:
                fields = line.rstrip('\n').split('\t')
                if len(fields) >= 6 and fields[3] == gene_name:
                    strand = fields[5]
                    break
        if strand is None:
            raise ValueError(
                f"Gene '{gene_name}' not found in {exons_bed} or {introns_bed}."
            )
        # No exon structure available; return empty list
        return None, None, None, strand, [], None

    exon_list = sorted(exons)
    gene_start = min(e[0] for e in exon_list)
    gene_end   = max(e[1] for e in exon_list)

    chrom_cm = NC_TO_CM.get(nc_chrom)
    if chrom_cm is None:
        raise ValueError(
            f"Cannot map chromosome '{nc_chrom}' to a CM accession. "
            "Add it to NC_TO_CM in the script."
        )

    return chrom_cm, gene_start, gene_end, strand, exon_list, nc_chrom


def load_intron_intervals(gene_name, introns_bed, chrom_cm):
    """
    Return a sorted list of (start, end) intron intervals for gene_name
    on chrom_cm, from a BED6 file using CM accessions.
    """
    intervals = []
    with open(introns_bed) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue
            if fields[3] != gene_name or fields[0] != chrom_cm:
                continue
            intervals.append((int(fields[1]), int(fields[2])))
    return sorted(intervals)


def load_mane_select_exons(gene_name, gtf_path, nc_chrom):
    """
    Return sorted list of (start, end) exon intervals from the MANE Select
    transcript for gene_name on nc_chrom, parsed from a GTF file.
    GTF is 1-based closed [start, end] → 0-based half-open [start-1, end).
    Returns None if no MANE Select transcript is found (caller should fall back).
    """
    exons = set()
    result = subprocess.run(
        ['grep', f'gene "{gene_name}"', gtf_path],
        capture_output=True, text=True
    )
    for line in result.stdout.splitlines():
        fields = line.split('\t')
        if len(fields) < 9:
            continue
        if fields[0] != nc_chrom or fields[2] != 'exon':
            continue
        if 'MANE Select' not in fields[8]:
            continue
        exons.add((int(fields[3]) - 1, int(fields[4])))
    return sorted(exons) if exons else None


def load_insertions(bwt_paths, chrom_cm, gene_start, gene_end, gene_strand,
                    intron_intervals):
    """
    Stream BWT files and collect intronic insertion positions within the gene.

    The .mapped_unique.bwt files are already deduplicated, so counts are
    summed directly across replicates. Only positions falling within
    intron_intervals are counted.

    Returns:
        sense_positions    — sorted list (same strand as gene)
        antisense_positions — sorted list (opposite strand from gene)
    """
    sense_all = []
    antisense_all = []

    for bwt_path in bwt_paths:
        if not os.path.exists(bwt_path):
            print(f"  [warn] BWT not found: {bwt_path}", file=sys.stderr)
            continue

        rep_label = os.path.basename(os.path.dirname(bwt_path))
        sense_rep = []
        antisense_rep = []

        with open(bwt_path) as fh:
            for line in fh:
                # BWT columns: read_name | strand | chrom | pos | seq | qual | n_other | mismatches
                fields = line.split('\t')
                if len(fields) < 5:
                    continue
                if fields[2] != chrom_cm:
                    continue
                bwt_strand = fields[1]
                bwt_pos    = int(fields[3])
                # Match the pipeline's BWT→BED coordinate: integration site is
                # at the LTR-genome junction (5' end of the read in genomic space).
                # For + strand reads this is bwt_pos; for - strand reads it is
                # bwt_pos + read_length - 1  (see analyze_sli.sh lines 316/323).
                if bwt_strand == '-':
                    site = bwt_pos + len(fields[4]) - 1
                else:
                    site = bwt_pos
                if site < gene_start or site > gene_end:
                    continue
                # Restrict to intronic positions only
                iv_idx = bisect.bisect_right(intron_intervals, (site, site)) - 1
                intronic = False
                for i in (iv_idx, iv_idx + 1):
                    if 0 <= i < len(intron_intervals):
                        s, e = intron_intervals[i]
                        if s <= site < e:
                            intronic = True
                            break
                if not intronic:
                    continue
                if bwt_strand == gene_strand:
                    sense_rep.append(site)
                else:
                    antisense_rep.append(site)

        print(
            f"  {rep_label}: sense={len(sense_rep):,}  antisense={len(antisense_rep):,}",
            file=sys.stderr
        )
        sense_all.extend(sense_rep)
        antisense_all.extend(antisense_rep)

    return sorted(sense_all), sorted(antisense_all)


# ══════════════════════════════════════════════════════════════════════════════
# Plotting helpers
# ══════════════════════════════════════════════════════════════════════════════

def bin_positions(positions, x_min, x_max, bin_size):
    """Return (bin_centers, counts) over [x_min, x_max]."""
    bins = np.arange(x_min, x_max + bin_size, bin_size)
    counts, edges = np.histogram(positions, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, counts


def draw_gene_diagram(ax, exon_list, x_left, x_right, gene_strand,
                      coord_fn=None):
    """
    Draw gene backbone + exon rectangles on ax.
    coord_fn(genomic_pos) → x-axis value; defaults to identity.
    x_left / x_right: limits of the x-axis in plot coordinates.
    """
    if coord_fn is None:
        coord_fn = lambda p: p

    ax.set_xlim(x_left, x_right)
    ax.set_ylim(0, 1)

    # Backbone
    ax.axhline(0.5, color='black', lw=1.5, zorder=1, solid_capstyle='butt')

    # Exons
    for estart, eend in exon_list:
        ex_l = min(coord_fn(estart), coord_fn(eend))
        ex_r = max(coord_fn(estart), coord_fn(eend))
        rect = mpatches.Rectangle(
            (ex_l, 0.15), ex_r - ex_l, 0.70,
            facecolor=EXON_FACE, edgecolor=EXON_EDGE, lw=0.5, zorder=2
        )
        ax.add_patch(rect)

    # Strand arrow at the 3' end of the gene (x_right is always 3' in plot coords)
    arrow_head = x_right
    arrow_tail = x_right + 0.04 * (x_left - x_right)
    ax.annotate(
        '', xy=(arrow_head, 0.5), xytext=(arrow_tail, 0.5),
        arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
    )

    ax.set_yticks([])
    ax.tick_params(bottom=False, top=False, left=False, right=False, labelbottom=False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _style_hist_axis(ax, top=True):
    """Remove tick marks and spines on the side adjacent to the gene diagram."""
    if top:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False)
    else:
        ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_locator(plt.MaxNLocator(4, integer=True))


# ══════════════════════════════════════════════════════════════════════════════
# Single-gene figure (Mode A)
# ══════════════════════════════════════════════════════════════════════════════

def plot_single_gene(gene_name, sense_pos, antisense_pos, diagram_exons,
                     gene_start, gene_end, gene_strand, bin_size):
    """
    Three-panel figure:
      Top    — sense insertions (bars grow upward away from gene)
      Middle — gene structure diagram
      Bottom — antisense insertions (bars grow downward away from gene)
    Always displayed in sense orientation (5' left, 3' right).
    Both histograms share the same y-scale.
    """
    centers, s_counts = bin_positions(sense_pos, gene_start, gene_end, bin_size)
    _,       a_counts = bin_positions(antisense_pos, gene_start, gene_end, bin_size)

    y_max = int(max(s_counts.max() if s_counts.size else 0,
                    a_counts.max() if a_counts.size else 0))

    # Sense orientation: for - strand genes put 5' end (high coord) on left
    x_plot_left  = gene_end   if gene_strand == '-' else gene_start
    x_plot_right = gene_start if gene_strand == '-' else gene_end

    fig = plt.figure(figsize=(18, 7))
    gs = GridSpec(3, 1, height_ratios=[2, 0.7, 2], hspace=0.05, figure=fig)

    ax_s = fig.add_subplot(gs[0])
    ax_g = fig.add_subplot(gs[1], sharex=ax_s)
    ax_a = fig.add_subplot(gs[2], sharex=ax_s)

    # ── Sense (top) ──────────────────────────────────────────────────────────
    ax_s.bar(centers, s_counts, width=bin_size * 0.9,
             color=SENSE_COLOR, alpha=0.85, linewidth=0)
    ax_s.set_ylim(0, y_max)
    ax_s.set_ylabel('Sense\ninsertions', fontsize=9)
    _style_hist_axis(ax_s, top=True)
    plt.setp(ax_s.get_xticklabels(), visible=False)
    ax_s.set_xlim(x_plot_left, x_plot_right)

    # ── Gene diagram (middle) ─────────────────────────────────────────────────
    draw_gene_diagram(ax_g, diagram_exons, x_plot_left, x_plot_right, gene_strand)
    plt.setp(ax_g.get_xticklabels(), visible=False)

    # ── Antisense (bottom) ────────────────────────────────────────────────────
    ax_a.bar(centers, a_counts, width=bin_size * 0.9,
             color=ANTISENSE_COLOR, alpha=0.85, linewidth=0)
    ax_a.set_ylim(0, y_max)
    ax_a.invert_yaxis()   # bars grow downward away from gene
    ax_a.set_ylabel('Antisense\ninsertions', fontsize=9)
    _style_hist_axis(ax_a, top=False)
    ax_a.set_xlabel('Genomic coordinate (bp)', fontsize=9)
    ax_a.set_xlim(x_plot_left, x_plot_right)
    ax_a.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x):,}'))

    # ── Title ─────────────────────────────────────────────────────────────────
    n_s = len(sense_pos)
    n_a = len(antisense_pos)
    fig.suptitle(
        f'{gene_name}  —  Retroviral insertion site distribution\n'
        f'Sense: {n_s:,} insertions   |   Antisense: {n_a:,} insertions   |   '
        f'Bin: {bin_size:,} bp   |   All replicates combined',
        fontsize=11, y=0.99
    )

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Multi-gene aligned figure (Modes B2 & B3)
# ══════════════════════════════════════════════════════════════════════════════

def _gene_to_plot_offset(pos, gene_start, gene_end, gene_strand):
    """Convert genomic position → gene-relative offset (0 = 5' end of gene)."""
    if gene_strand == '-':
        return gene_end - pos
    return pos - gene_start


def plot_aligned(gene_data_list, alignment, bin_size, title):
    """
    Vertically stacked panels, one per gene, sharing a common x-axis.

    alignment: '5prime'          → x=0 is the 5' end (TSS)
               'antisense_peak'  → x=0 is the peak antisense bin
    """
    n = len(gene_data_list)

    # ── First pass: compute per-gene offsets and global x range ──────────────
    ref_offsets = []
    all_x_mins, all_x_maxs = [], []

    for gd in gene_data_list:
        gs, ge, gst = gd['gene_start'], gd['gene_end'], gd['gene_strand']
        gene_len = ge - gs

        anti_gene = [_gene_to_plot_offset(p, gs, ge, gst) for p in gd['antisense_pos']]
        gd['_sense_gene']  = [_gene_to_plot_offset(p, gs, ge, gst) for p in gd['sense_pos']]
        gd['_anti_gene']   = anti_gene
        gd['_gene_len']    = gene_len

        if alignment == '5prime':
            ref = 0
        else:
            if anti_gene:
                bins = np.arange(0, gene_len + bin_size, bin_size)
                counts, edges = np.histogram(anti_gene, bins=bins)
                pk = np.argmax(counts)
                ref = (edges[pk] + edges[pk + 1]) / 2
            else:
                ref = 0

        ref_offsets.append(ref)
        all_x_mins.append(0 - ref)
        all_x_maxs.append(gene_len - ref)

    x_min = min(all_x_mins) - bin_size
    x_max = max(all_x_maxs) + bin_size

    # ── Build figure ──────────────────────────────────────────────────────────
    row_heights = [1.8, 0.5, 1.8] * n
    fig_h = max(4 * n, 8)
    fig = plt.figure(figsize=(18, fig_h))
    gs_layout = GridSpec(n * 3, 1, height_ratios=row_heights, hspace=0.05, figure=fig)

    first_ax = None

    for i, gd in enumerate(gene_data_list):
        ref = ref_offsets[i]

        # Shift positions into alignment-relative coordinates
        sense_x = [o - ref for o in gd['_sense_gene']]
        anti_x  = [o - ref for o in gd['_anti_gene']]

        centers, s_counts = bin_positions(sense_x, x_min, x_max, bin_size)
        _,       a_counts = bin_positions(anti_x,  x_min, x_max, bin_size)

        base = i * 3
        if first_ax is None:
            ax_s = fig.add_subplot(gs_layout[base])
            ax_g = fig.add_subplot(gs_layout[base + 1], sharex=ax_s)
            ax_a = fig.add_subplot(gs_layout[base + 2], sharex=ax_s)
            first_ax = ax_s
        else:
            ax_s = fig.add_subplot(gs_layout[base],     sharex=first_ax)
            ax_g = fig.add_subplot(gs_layout[base + 1], sharex=first_ax)
            ax_a = fig.add_subplot(gs_layout[base + 2], sharex=first_ax)

        # Sense
        y_max = int(max(s_counts.max() if s_counts.size else 0,
                        a_counts.max() if a_counts.size else 0))

        ax_s.bar(centers, s_counts, width=bin_size * 0.9,
                 color=SENSE_COLOR, alpha=0.85, linewidth=0)
        ax_s.set_ylim(0, y_max)
        ax_s.set_ylabel(gd['gene_name'], fontsize=9, rotation=0,
                        labelpad=50, va='center', ha='right')
        ax_s.tick_params(labelbottom=False, bottom=False)
        ax_s.spines['bottom'].set_visible(False)
        ax_s.spines['right'].set_visible(False)
        ax_s.yaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
        ax_s.set_xlim(x_min, x_max)

        # Gene diagram (in alignment-relative coords)
        gene_x_start = 0 - ref
        gene_x_end   = gd['_gene_len'] - ref
        x_l = min(gene_x_start, gene_x_end)
        x_r = max(gene_x_start, gene_x_end)

        gs_, ge_, gst_ = gd['gene_start'], gd['gene_end'], gd['gene_strand']

        def make_coord_fn(gs_=gs_, ge_=ge_, gst_=gst_, ref_=ref):
            def coord_fn(p):
                return _gene_to_plot_offset(p, gs_, ge_, gst_) - ref_
            return coord_fn

        draw_gene_diagram(ax_g, gd['diagram_exons'], x_min, x_max, gst_,
                          coord_fn=make_coord_fn())
        ax_g.tick_params(labelbottom=False, bottom=False)

        # Antisense
        ax_a.bar(centers, a_counts, width=bin_size * 0.9,
                 color=ANTISENSE_COLOR, alpha=0.85, linewidth=0)
        ax_a.set_ylim(0, y_max)
        ax_a.invert_yaxis()   # bars grow downward away from gene
        ax_a.spines['top'].set_visible(False)
        ax_a.spines['right'].set_visible(False)
        ax_a.yaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
        ax_a.set_xlim(x_min, x_max)

        if i < n - 1:
            ax_a.tick_params(labelbottom=False)
        else:
            if alignment == '5prime':
                ax_a.set_xlabel("Offset from 5' end (bp)", fontsize=9)
            else:
                ax_a.set_xlabel("Offset from peak antisense bin (bp)", fontsize=9)

    # Legend
    legend_handles = [
        mpatches.Patch(facecolor=SENSE_COLOR,     alpha=0.85, label='Sense insertions'),
        mpatches.Patch(facecolor=ANTISENSE_COLOR,  alpha=0.85, label='Antisense insertions'),
        mpatches.Patch(facecolor=EXON_FACE,        alpha=1.0,  label='Exon'),
    ]
    fig.legend(handles=legend_handles, loc='upper right', fontsize=8,
               bbox_to_anchor=(1.0, 1.0))

    fig.suptitle(title, fontsize=11, y=1.01)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--gene', metavar='NAME',
                      help='Single gene name (e.g. STX18)')
    mode.add_argument('--input', metavar='CSV',
                      help='sense_vs_antisense_counts.csv (must have gene and p columns)')

    parser.add_argument('--top', type=int, default=10, metavar='N',
                        help='Number of top genes by p-value to use (default: 10). '
                             'Only applies with --input.')
    parser.add_argument('--bwt', nargs='+', metavar='PATH',
                        help='BWT alignment file(s). Accepts shell glob patterns. '
                             'Default: ControlData replicate 1–4.')
    parser.add_argument('--exons',   default=DEFAULT_EXONS,   metavar='PATH',
                        help='Exons BED file (NC accessions).')
    parser.add_argument('--introns', default=DEFAULT_INTRONS, metavar='PATH',
                        help='Introns BED file (CM accessions, fallback for gene lookup).')
    parser.add_argument('--gtf', default=DEFAULT_GTF, metavar='PATH',
                        help='GTF file for MANE Select canonical exon diagram.')
    parser.add_argument('--bin-size', type=int, default=500, metavar='BP',
                        help='Histogram bin width in base pairs (default: 500).')
    parser.add_argument('--output-dir', default='/scratch/alpine/wist9668/figures',
                        metavar='DIR',
                        help='Directory for output PDFs '
                             '(default: /scratch/alpine/wist9668/figures).')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Resolve BWT paths ─────────────────────────────────────────────────────
    if args.bwt:
        bwt_paths = []
        for pattern in args.bwt:
            expanded = sorted(glob.glob(pattern))
            bwt_paths.extend(expanded if expanded else [pattern])
    else:
        bwt_paths = DEFAULT_BWT_PATHS

    # ── Determine gene list ───────────────────────────────────────────────────
    if args.gene:
        gene_list = [args.gene]
    else:
        rows = []
        with open(args.input, newline='') as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                try:
                    rows.append((row['gene'], float(row['p'])))
                except (KeyError, ValueError):
                    continue
        rows.sort(key=lambda x: x[1])
        gene_list = [g for g, _ in rows[:args.top]]
        print(f"Top {args.top} genes by p-value: {', '.join(gene_list)}", file=sys.stderr)

    # ── Load gene structures and insertions ───────────────────────────────────
    gene_data_list = []

    for gene_name in gene_list:
        print(f"\n{'─'*50}", file=sys.stderr)
        print(f"Processing {gene_name} ...", file=sys.stderr)

        try:
            chrom_cm, g_start, g_end, g_strand, exon_list, nc_chrom = load_gene_structure(
                gene_name, args.exons, args.introns
            )
        except ValueError as err:
            print(f"  Skipping: {err}", file=sys.stderr)
            continue

        if chrom_cm is None:
            print(f"  Skipping {gene_name}: no chromosome mapping available.", file=sys.stderr)
            continue

        intron_ivs = load_intron_intervals(gene_name, args.introns, chrom_cm)
        if not intron_ivs:
            print(f"  Skipping {gene_name}: no intron intervals found in {args.introns}.",
                  file=sys.stderr)
            continue

        mane_exons = load_mane_select_exons(gene_name, args.gtf, nc_chrom)
        diagram_exons = mane_exons if mane_exons else exon_list
        mane_label = f"MANE Select ({len(mane_exons)} exons)" if mane_exons else "all isoforms"

        print(
            f"  Location : {chrom_cm}:{g_start:,}–{g_end:,} ({g_strand} strand)\n"
            f"  Gene span: {g_end - g_start:,} bp   Diagram: {mane_label}",
            file=sys.stderr
        )

        sense_pos, anti_pos = load_insertions(
            bwt_paths, chrom_cm, g_start, g_end, g_strand, intron_ivs
        )
        print(
            f"  Total    : sense={len(sense_pos):,}  antisense={len(anti_pos):,}",
            file=sys.stderr
        )

        gene_data_list.append({
            'gene_name':     gene_name,
            'sense_pos':     sense_pos,
            'antisense_pos': anti_pos,
            'exon_list':     exon_list,
            'diagram_exons': diagram_exons,
            'gene_start':    g_start,
            'gene_end':      g_end,
            'gene_strand':   g_strand,
        })

    if not gene_data_list:
        print("\nError: no valid genes to plot.", file=sys.stderr)
        sys.exit(1)

    bin_size = args.bin_size
    print(f"\n{'─'*50}", file=sys.stderr)

    # ── Single-gene mode ──────────────────────────────────────────────────────
    if args.gene:
        gd  = gene_data_list[0]
        fig = plot_single_gene(
            gd['gene_name'], gd['sense_pos'], gd['antisense_pos'],
            gd['diagram_exons'], gd['gene_start'], gd['gene_end'], gd['gene_strand'],
            bin_size
        )
        out = os.path.join(args.output_dir, f"{gd['gene_name']}_insertion_map.pdf")
        fig.savefig(out, bbox_inches='tight', dpi=150)
        plt.close(fig)
        print(f"Saved: {out}", file=sys.stderr)

    # ── Multi-gene mode ───────────────────────────────────────────────────────
    else:
        screen = os.path.splitext(os.path.basename(args.input))[0]
        top_n  = len(gene_data_list)

        # B1: Individual plots in one multi-page PDF
        out_indiv = os.path.join(args.output_dir,
                                 f"{screen}_top{top_n}_individual.pdf")
        with PdfPages(out_indiv) as pdf:
            for gd in gene_data_list:
                fig = plot_single_gene(
                    gd['gene_name'], gd['sense_pos'], gd['antisense_pos'],
                    gd['diagram_exons'], gd['gene_start'], gd['gene_end'],
                    gd['gene_strand'], bin_size
                )
                pdf.savefig(fig, bbox_inches='tight', dpi=150)
                plt.close(fig)
        print(f"Saved: {out_indiv}", file=sys.stderr)

        # B2: Stacked panels aligned by 5' end
        out_5p = os.path.join(args.output_dir,
                              f"{screen}_top{top_n}_aligned_5prime.pdf")
        fig = plot_aligned(
            gene_data_list,
            alignment='5prime',
            bin_size=bin_size,
            title=(f"Top {top_n} genes — aligned by 5' end   "
                   f"(bin size: {bin_size:,} bp, all replicates combined)")
        )
        fig.savefig(out_5p, bbox_inches='tight', dpi=150)
        plt.close(fig)
        print(f"Saved: {out_5p}", file=sys.stderr)

        # B3: Stacked panels aligned by peak antisense bin
        out_ap = os.path.join(args.output_dir,
                              f"{screen}_top{top_n}_aligned_antisense_peak.pdf")
        fig = plot_aligned(
            gene_data_list,
            alignment='antisense_peak',
            bin_size=bin_size,
            title=(f"Top {top_n} genes — aligned by peak antisense bin   "
                   f"(bin size: {bin_size:,} bp, all replicates combined)")
        )
        fig.savefig(out_ap, bbox_inches='tight', dpi=150)
        plt.close(fig)
        print(f"Saved: {out_ap}", file=sys.stderr)


if __name__ == '__main__':
    main()
