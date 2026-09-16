#!/usr/bin/env python3
"""
Count unique exonic insertion sites per gene from a bowtie BWT file.

Exonic insertions on either strand are inactivating (gene-trap or frameshift),
so strand is ignored when assigning insertions to genes — but (chrom, strand, pos)
is still used as the deduplication key to match the existing pipeline convention.

Uses bedtools intersect internally. Requires bedtools on PATH.

Usage:
    python3 count_exon_insertions.py --bwt <file> --exons <file> --out <file>

Output: TSV with columns: gene, exon_insertions
"""
import argparse
import subprocess
import tempfile
import os
import sys
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--bwt',   required=True, help='replicate.mapped_unique.bwt')
    p.add_argument('--exons', required=True, help='exons_CM_mapped.bed (6-col BED)')
    p.add_argument('--out',   required=True, help='Output TSV path')
    return p.parse_args()

def bwt_to_bed_dedup(bwt_path):
    """Read BWT, deduplicate by (chrom, strand, pos), return sorted BED6 lines."""
    seen = set()
    lines = []
    with open(bwt_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            strand = parts[1]
            chrom  = parts[2]
            pos    = parts[3]
            key = (chrom, strand, pos)
            if key not in seen:
                seen.add(key)
                lines.append(f"{chrom}\t{pos}\t{int(pos)+1}\t.\t.\t{strand}\n")
    lines.sort()
    return lines

def main():
    args = parse_args()

    # Write deduplicated insertion sites to a temp BED file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
        tmp_path = tmp.name
        for line in bwt_to_bed_dedup(args.bwt):
            tmp.write(line)

    try:
        # bedtools intersect: for each insertion site, report overlapping exon entry
        # exons BED col 3 = gene name  (0-indexed)
        # -wb reports both the a-side (insertion) and b-side (exon) columns
        result = subprocess.run(
            ['bedtools', 'intersect', '-a', tmp_path, '-b', args.exons, '-wb'],
            capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"bedtools error: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    finally:
        os.unlink(tmp_path)

    # Parse bedtools output:
    # cols 0-5: insertion BED6 (chrom, start, end, ., ., strand)
    # cols 6-11: exon BED6 (chrom, start, end, gene_name, score, gene_strand)
    # gene_name is col 9
    counts = defaultdict(set)
    for line in result.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) < 10:
            continue
        chrom, start, _, _, _, strand = parts[:6]
        gene = parts[9]
        counts[gene].add((chrom, start, strand))

    with open(args.out, 'w') as f:
        f.write("gene\texon_insertions\n")
        for gene, sites in sorted(counts.items()):
            f.write(f"{gene}\t{len(sites)}\n")

    n_genes = len(counts)
    n_sites = sum(len(s) for s in counts.values())
    print(f"Done: {n_sites} unique exonic insertion sites across {n_genes} genes",
          file=sys.stderr)

if __name__ == '__main__':
    main()
