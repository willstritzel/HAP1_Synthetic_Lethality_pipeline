#!/usr/bin/env python3
"""
detect_strand_bias.py

Sliding window binomial test for strand bias across the genome.

Input:  BED6 from aggregate_insertions.py (CM chromosome coordinates)
Output: 9-column BED of significant strand-biased regions after genome-wide
        Benjamini-Hochberg FDR correction and merging of adjacent windows.

Output columns (tab-separated, header line starts with '#'):
  chrom, start, end, direction, sense_count, antisense_count,
  sense_fraction, p_value, fdr

  direction: sense_bias (+ strand enriched) or antisense_bias (- strand enriched)
"""

import argparse
import os
import sys

import numpy as np
import scipy.stats
import statsmodels.sandbox.stats.multicomp

sys.path.insert(0, os.path.dirname(__file__))
from make_tiled_windows import CHROM_ORDER, GRCH38_SIZES  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, help="BED6 from aggregate_insertions.py")
    p.add_argument("--output", required=True, help="Output biased regions BED")
    p.add_argument(
        "--window-size",
        type=int,
        default=10000,
        help="Sliding window size in bp (default: 10000)",
    )
    p.add_argument(
        "--step-size",
        type=int,
        default=1000,
        help="Sliding window step in bp (default: 1000)",
    )
    p.add_argument(
        "--min-insertions",
        type=int,
        default=10,
        help="Min total insertions in a window to run the binomial test (default: 10)",
    )
    p.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="BH FDR significance threshold (default: 0.05)",
    )
    p.add_argument(
        "--merge-gap",
        type=int,
        default=1000,
        help="Max bp gap between adjacent significant windows to merge (default: 1000)",
    )
    return p.parse_args()


def load_insertions(bed_path, chrom_order):
    """
    Read BED6 and build per-chromosome sorted numpy int64 position arrays.
    Only primary chromosomes (in chrom_order) are loaded; alt contigs are skipped.
    Returns (sense_pos, antisense_pos) where each is dict[chrom -> np.ndarray].
    """
    chrom_set = set(chrom_order)
    sense_lists = {c: [] for c in chrom_order}
    antisense_lists = {c: [] for c in chrom_order}
    skipped = set()

    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom = parts[0]
            if chrom not in chrom_set:
                skipped.add(chrom)
                continue
            pos = int(parts[1])
            strand = parts[5].strip()
            if strand == "+":
                sense_lists[chrom].append(pos)
            else:
                antisense_lists[chrom].append(pos)

    if skipped:
        print(f"  Skipped {len(skipped)} alt/unplaced contig(s)", file=sys.stderr)

    sense_pos = {
        c: np.array(sorted(sense_lists[c]), dtype=np.int64) for c in chrom_order
    }
    antisense_pos = {
        c: np.array(sorted(antisense_lists[c]), dtype=np.int64) for c in chrom_order
    }
    return sense_pos, antisense_pos


def count_in_window(sorted_arr, start, end):
    """Count elements in sorted_arr within [start, end) using binary search."""
    lo = np.searchsorted(sorted_arr, start, side="left")
    hi = np.searchsorted(sorted_arr, end, side="left")
    return int(hi - lo)


def scan_chromosome(chrom, chrom_size, sense_arr, antisense_arr,
                    window_size, step_size, min_insertions):
    """
    Slide a window across one chromosome and collect per-window stats.

    Returns list of (chrom, start, end, sense_count, antisense_count, p_value).
    Windows where total < min_insertions get p_value=None (excluded from FDR).
    """
    results = []
    start = 0
    while start < chrom_size:
        end = min(start + window_size, chrom_size)
        sense = count_in_window(sense_arr, start, end)
        antisense = count_in_window(antisense_arr, start, end)
        total = sense + antisense
        if total >= min_insertions:
            p = scipy.stats.binomtest(
                sense, total, p=0.5, alternative="two-sided"
            ).pvalue
            results.append((chrom, start, end, sense, antisense, p))
        else:
            results.append((chrom, start, end, sense, antisense, None))
        start += step_size
    return results


def apply_fdr_correction(all_windows, fdr_threshold):
    """
    Genome-wide BH FDR correction across all testable windows.

    Collects all non-None p-values into one list, applies multipletests()
    (same pattern as binomial_test.py), then returns only windows with FDR
    <= fdr_threshold, with direction label added.

    Returns list of (chrom, start, end, sense, antisense, p_value, fdr, direction).
    """
    testable_idx = [i for i, w in enumerate(all_windows) if w[5] is not None]
    pvalues = [all_windows[i][5] for i in testable_idx]

    if not pvalues:
        print(
            "  WARNING: No windows met the minimum insertion threshold.",
            file=sys.stderr,
        )
        return []

    _, fdr_values, _, _ = statsmodels.sandbox.stats.multicomp.multipletests(
        pvalues, alpha=fdr_threshold, method="fdr_bh", returnsorted=False
    )

    sig = []
    for rank, idx in enumerate(testable_idx):
        fdr = fdr_values[rank]
        if fdr <= fdr_threshold:
            chrom, start, end, sense, antisense, p = all_windows[idx]
            direction = "sense_bias" if sense > antisense else "antisense_bias"
            sig.append((chrom, start, end, sense, antisense, p, fdr, direction))

    return sig


def merge_significant_windows(sig_windows, merge_gap, chrom_order):
    """
    Merge adjacent significant windows sharing the same chromosome and direction.

    Two windows are merged when the gap between them <= merge_gap AND they have
    the same direction. Merged regions sum sense/antisense counts and take the
    minimum p-value and FDR of constituent windows.

    Returns sorted list of merged regions.
    """
    if not sig_windows:
        return []

    chrom_idx = {c: i for i, c in enumerate(chrom_order)}
    sig_windows = sorted(
        sig_windows, key=lambda w: (chrom_idx.get(w[0], 9999), w[1])
    )

    merged = []
    cur = list(sig_windows[0])

    for w in sig_windows[1:]:
        chrom, start, end, sense, antisense, p, fdr, direction = w
        same_chrom = chrom == cur[0]
        same_dir = direction == cur[7]
        adjacent = start - cur[2] <= merge_gap

        if same_chrom and same_dir and adjacent:
            cur[2] = end
            cur[3] += sense
            cur[4] += antisense
            cur[5] = min(cur[5], p)
            cur[6] = min(cur[6], fdr)
        else:
            merged.append(tuple(cur))
            cur = list(w)

    merged.append(tuple(cur))
    return merged


def write_output(merged_regions, outpath):
    """Write 9-column BED with a comment header line."""
    with open(outpath, "w") as out:
        out.write(
            "#chrom\tstart\tend\tdirection\tsense_count\tantisense_count\t"
            "sense_fraction\tp_value\tfdr\n"
        )
        for chrom, start, end, sense, antisense, p, fdr, direction in merged_regions:
            total = sense + antisense
            frac = sense / total if total > 0 else 0.0
            out.write(
                f"{chrom}\t{start}\t{end}\t{direction}\t{sense}\t{antisense}\t"
                f"{frac:.4f}\t{p:.6e}\t{fdr:.6e}\n"
            )


def main():
    args = parse_args()

    print(f"[detect_strand_bias] Input:            {args.input}", file=sys.stderr)
    print(f"[detect_strand_bias] Output:           {args.output}", file=sys.stderr)
    print(f"[detect_strand_bias] Window size:      {args.window_size:,} bp", file=sys.stderr)
    print(f"[detect_strand_bias] Step size:        {args.step_size:,} bp", file=sys.stderr)
    print(f"[detect_strand_bias] Min insertions:   {args.min_insertions}", file=sys.stderr)
    print(f"[detect_strand_bias] FDR threshold:    {args.fdr_threshold}", file=sys.stderr)
    print(f"[detect_strand_bias] Merge gap:        {args.merge_gap:,} bp\n", file=sys.stderr)

    print("[detect_strand_bias] Loading insertions...", file=sys.stderr)
    sense_pos, antisense_pos = load_insertions(args.input, CHROM_ORDER)

    all_windows = []
    for chrom in CHROM_ORDER:
        chrom_size = GRCH38_SIZES[chrom]
        s_arr = sense_pos[chrom]
        a_arr = antisense_pos[chrom]
        print(
            f"  {chrom}: {len(s_arr):,} sense, {len(a_arr):,} antisense insertions",
            file=sys.stderr,
        )
        windows = scan_chromosome(
            chrom, chrom_size, s_arr, a_arr,
            args.window_size, args.step_size, args.min_insertions,
        )
        all_windows.extend(windows)

    testable = sum(1 for w in all_windows if w[5] is not None)
    print(
        f"\n[detect_strand_bias] Total windows:    {len(all_windows):,}",
        file=sys.stderr,
    )
    print(
        f"[detect_strand_bias] Testable windows: {testable:,}", file=sys.stderr
    )

    print(
        "\n[detect_strand_bias] Applying genome-wide BH FDR correction...",
        file=sys.stderr,
    )
    sig_windows = apply_fdr_correction(all_windows, args.fdr_threshold)
    print(f"  Significant windows: {len(sig_windows):,}", file=sys.stderr)

    print(
        f"[detect_strand_bias] Merging adjacent windows (gap <= {args.merge_gap} bp)...",
        file=sys.stderr,
    )
    merged = merge_significant_windows(sig_windows, args.merge_gap, CHROM_ORDER)
    sense_biased = sum(1 for r in merged if r[7] == "sense_bias")
    antisense_biased = len(merged) - sense_biased
    print(f"  Merged regions:      {len(merged):,}", file=sys.stderr)
    print(f"  Sense-biased:        {sense_biased:,}", file=sys.stderr)
    print(f"  Antisense-biased:    {antisense_biased:,}", file=sys.stderr)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    write_output(merged, args.output)
    print(f"\n[detect_strand_bias] Done. Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
