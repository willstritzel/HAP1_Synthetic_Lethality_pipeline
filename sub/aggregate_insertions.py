#!/usr/bin/env python3
"""
aggregate_insertions.py

Pools raw bowtie alignment files from all 4 replicates of a tiled HAP1 screen,
de-duplicates insertion sites across replicates, and outputs a sorted BED6 file.

Each line in the output represents one unique genomic insertion site.

BED6 columns: chrom, start, end, bitmask, n_replicates, strand
  - chrom/start/end: integration site in 0-based half-open coordinates (CM accession space)
  - bitmask:         integer where bit i is set if replicate (i+1) observed this site
  - n_replicates:    popcount(bitmask), i.e. 1-4
  - strand:          + or -

The dedup key (chrom, strand, pos) matches the existing pipeline convention
used by analyze_sli.sh: awk '!seen[$3$2$4]++'
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from make_tiled_windows import CHROM_ORDER, GRCH38_SIZES  # noqa: E402  (after sys.path)

DEFAULT_REP_DIR = (
    "/scratch/alpine/wist9668/data/analyzed_data/"
    "synthetic_lethal_screens/tiled_2kb_hap1-4_output"
)


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--rep-dir",
        default=DEFAULT_REP_DIR,
        help="Directory containing replicate_{1..4}/ subdirectories (default: %(default)s)",
    )
    p.add_argument("--output", required=True, help="Output BED6 file path")
    p.add_argument(
        "--screen-name",
        default="hap1-4_tiled",
        help="Screen name for progress messages (default: %(default)s)",
    )
    return p.parse_args()


def discover_bwt_files(rep_dir):
    """Return [(rep_number, filepath), ...] for all found replicate .bwt files."""
    found = []
    for rep_num in range(1, 5):
        path = os.path.join(rep_dir, f"replicate_{rep_num}", "replicate.mapped_unique.bwt")
        if os.path.isfile(path):
            found.append((rep_num, path))
        else:
            print(f"  WARNING: Not found: {path}", file=sys.stderr)
    if not found:
        raise FileNotFoundError(
            f"No replicate.mapped_unique.bwt files found under {rep_dir}"
        )
    return found


def load_bwt_file(filepath, rep_idx, insertions):
    """
    Stream one bowtie output file and update insertions dict in place.

    BWT columns (tab-separated):
      0: read_name  1: strand  2: chrom  3: pos (0-based)  4: seq  5: qual  ...

    Sets bit (1 << rep_idx) in insertions[(chrom, strand, pos)].
    Returns the number of lines processed.
    """
    bit = 1 << rep_idx
    count = 0
    with open(filepath) as fh:
        for line in fh:
            parts = line.split("\t")
            chrom = sys.intern(parts[2])
            strand = parts[1]
            pos = int(parts[3])
            key = (chrom, strand, pos)
            if key in insertions:
                insertions[key] |= bit
            else:
                insertions[key] = bit
            count += 1
    return count


def write_bed6(insertions, outpath, chrom_order):
    """
    Write sorted BED6 to outpath.

    Sort order: canonical chromosomes first (CHROM_ORDER), then alt contigs
    alphabetically. Within each chromosome, sort by position ascending.

    Returns summary dict: {chrom: {'sense': int, 'antisense': int}}.
    """
    chrom_set = set(chrom_order)

    # Bucket all keys by chromosome
    by_chrom = {}
    for (chrom, strand, pos), bitmask in insertions.items():
        if chrom not in by_chrom:
            by_chrom[chrom] = []
        by_chrom[chrom].append((pos, strand, bitmask))

    canonical = [c for c in chrom_order if c in by_chrom]
    alt = sorted(c for c in by_chrom if c not in chrom_set)
    ordered = canonical + alt

    summary = {}
    with open(outpath, "w") as out:
        for chrom in ordered:
            entries = sorted(by_chrom[chrom])  # sort by pos (first tuple element)
            sense = antisense = 0
            for pos, strand, bitmask in entries:
                n_reps = bin(bitmask).count("1")
                out.write(f"{chrom}\t{pos}\t{pos + 1}\t{bitmask}\t{n_reps}\t{strand}\n")
                if strand == "+":
                    sense += 1
                else:
                    antisense += 1
            summary[chrom] = {"sense": sense, "antisense": antisense}

    return summary


def main():
    args = parse_args()

    print(f"[aggregate_insertions] Screen:   {args.screen_name}", file=sys.stderr)
    print(f"[aggregate_insertions] Rep dir:  {args.rep_dir}", file=sys.stderr)
    print(f"[aggregate_insertions] Output:   {args.output}", file=sys.stderr)

    bwt_files = discover_bwt_files(args.rep_dir)
    print(
        f"[aggregate_insertions] Found {len(bwt_files)} replicate(s)\n", file=sys.stderr
    )

    insertions = {}
    for rep_num, filepath in bwt_files:
        print(
            f"[aggregate_insertions] Loading replicate {rep_num}: {os.path.basename(filepath)}",
            file=sys.stderr,
        )
        n = load_bwt_file(filepath, rep_num - 1, insertions)
        print(f"  {n:,} lines processed", file=sys.stderr)

    total_unique = len(insertions)
    print(
        f"\n[aggregate_insertions] Unique sites after cross-replicate de-dup: {total_unique:,}",
        file=sys.stderr,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    print(f"[aggregate_insertions] Writing output...", file=sys.stderr)
    summary = write_bed6(insertions, args.output, CHROM_ORDER)

    # Per-chromosome summary table
    print(
        f"\n{'Chrom':<15} {'Sense (+ strand)':>18} {'Antisense (- strand)':>21} {'Total':>12}",
        file=sys.stderr,
    )
    print("-" * 70, file=sys.stderr)
    for chrom in CHROM_ORDER:
        if chrom in summary:
            s = summary[chrom]["sense"]
            a = summary[chrom]["antisense"]
            print(f"{chrom:<15} {s:>18,} {a:>21,} {s + a:>12,}", file=sys.stderr)

    print(f"\n[aggregate_insertions] Done. Total unique sites: {total_unique:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
