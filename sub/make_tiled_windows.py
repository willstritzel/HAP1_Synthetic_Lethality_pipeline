#!/usr/bin/env python3
"""
Generate tiled genomic window .bed files for annotation-independent haploid analysis.

Two bed files are created:
  - genome_10kb_windows.bed
        Non-overlapping 10kb windows starting at position 0.
        Window borders at: 0, 10000, 20000, 30000, ...

  - genome_2kb_windows_1kb_offset.bed
        Non-overlapping 2kb windows with a 1kb leading window to offset borders.
        First window: [0, 1000). All subsequent windows: 2000bp.
        Window borders at: 0, 1000, 3000, 5000, 7000, 9000, 11000, ...
        This ensures every 10kb border (10000, 20000, ...) falls exactly in
        the center of a 2kb window, 1000bp from either border — preventing
        border-alignment artifacts between the two bed files.

Strand is set to "+" for all windows so the pipeline correctly classifies
sense (+) and antisense (-) strand insertions within each window.

Output format (6 columns, tab-separated):
  chrom  start  end  window_name  .  +
"""

import os
import argparse

# GRCh38 primary assembly chromosome sizes (Ensembl release 109 / UCSC hg38)
# Chromosome names use GenBank/Ensembl accessions to match the bowtie reference index
# (GCA_000001405.15_GRCh38_genomic_bowtie) and the functional annotation BED files.
GRCH38_SIZES = {
    "CM000663.2": 248956422,  # chr1
    "CM000664.2": 242193529,  # chr2
    "CM000665.2": 198295559,  # chr3
    "CM000666.2": 190214555,  # chr4
    "CM000667.2": 181538259,  # chr5
    "CM000668.2": 170805979,  # chr6
    "CM000669.2": 159345973,  # chr7
    "CM000670.2": 145138636,  # chr8
    "CM000671.2": 138394717,  # chr9
    "CM000672.2": 133797422,  # chr10
    "CM000673.2": 135086622,  # chr11
    "CM000674.2": 133275309,  # chr12
    "CM000675.2": 114364328,  # chr13
    "CM000676.2": 107043718,  # chr14
    "CM000677.2": 101991189,  # chr15
    "CM000678.2":  90338345,  # chr16
    "CM000679.2":  83257441,  # chr17
    "CM000680.2":  80373285,  # chr18
    "CM000681.2":  58617616,  # chr19
    "CM000682.2":  64444167,  # chr20
    "CM000683.2":  46709983,  # chr21
    "CM000684.2":  50818468,  # chr22
    "CM000685.2": 156040895,  # chrX
    "CM000686.2":  57227415,  # chrY
}

CHROM_ORDER = [
    "CM000663.2", "CM000664.2", "CM000665.2", "CM000666.2", "CM000667.2", "CM000668.2",
    "CM000669.2", "CM000670.2", "CM000671.2", "CM000672.2", "CM000673.2", "CM000674.2",
    "CM000675.2", "CM000676.2", "CM000677.2", "CM000678.2", "CM000679.2", "CM000680.2",
    "CM000681.2", "CM000682.2", "CM000683.2", "CM000684.2", "CM000685.2", "CM000686.2",
]


def make_uniform_windows(chrom, chrom_size, window_size, outfile):
    """Write uniform non-overlapping windows for one chromosome."""
    count = 0
    start = 0
    while start < chrom_size:
        end = min(start + window_size, chrom_size)
        name = f"{chrom}:{start}-{end}"
        outfile.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t+\n")
        start += window_size
        count += 1
    return count


def make_offset_windows(chrom, chrom_size, window_size, lead_size, outfile):
    """Write offset windows: first window is lead_size, rest are window_size.

    This places uniform-window borders in the interior of these windows.
    E.g. with window_size=2000, lead_size=1000:
      [0, 1000), [1000, 3000), [3000, 5000), ...
      Borders at 1000, 3000, 5000, 7000, 9000, 11000, ...
      Every 10kb border (10000, 20000, ...) falls at the midpoint of a window.
    """
    count = 0
    # First window (lead)
    if lead_size > 0 and lead_size < chrom_size:
        end = lead_size
        name = f"{chrom}:0-{end}"
        outfile.write(f"{chrom}\t0\t{end}\t{name}\t.\t+\n")
        count += 1
        start = lead_size
    else:
        start = 0

    while start < chrom_size:
        end = min(start + window_size, chrom_size)
        name = f"{chrom}:{start}-{end}"
        outfile.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t+\n")
        start += window_size
        count += 1
    return count


def generate_bed(output_path, filename, window_size, lead_size=0):
    filepath = os.path.join(output_path, filename)
    total = 0
    with open(filepath, "w") as f:
        for chrom in CHROM_ORDER:
            chrom_size = GRCH38_SIZES[chrom]
            if lead_size:
                n = make_offset_windows(chrom, chrom_size, window_size, lead_size, f)
            else:
                n = make_uniform_windows(chrom, chrom_size, window_size, f)
            total += n
    print(f"Written {total} windows to {filepath}")
    return filepath


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        default="/scratch/alpine/wist9668/reference/annotation/tiled_windows",
        help="Output directory for .bed files",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # bed_A: uniform 10kb windows
    generate_bed(args.outdir,
                 filename="genome_10kb_windows.bed",
                 window_size=10000)

    # bed_B: 2kb windows with 1kb lead window to offset borders by 1kb
    generate_bed(args.outdir,
                 filename="genome_2kb_windows_1kb_offset.bed",
                 window_size=2000,
                 lead_size=1000)


if __name__ == "__main__":
    main()
