#!/usr/bin/env python3
"""
Filter lncRNA introns BED to remove entries for lncRNA genes that overlap
protein-coding gene bodies (sense or antisense).

Usage:
    python3 filter_lncRNA_introns.py <overlap_gene_list> <input_bed> <output_bed>
"""
import sys

if len(sys.argv) != 4:
    print("Usage: filter_lncRNA_introns.py <overlap_gene_list> <input_bed> <output_bed>",
          file=sys.stderr)
    sys.exit(1)

overlap_file, input_bed, output_bed = sys.argv[1], sys.argv[2], sys.argv[3]

with open(overlap_file) as f:
    exclude = set(line.strip() for line in f if line.strip())

kept = excluded = 0
with open(input_bed) as fin, open(output_bed, "w") as fout:
    for line in fin:
        gene = line.split("\t")[3]
        if gene not in exclude:
            fout.write(line)
            kept += 1
        else:
            excluded += 1

print(f"Kept: {kept}  Excluded: {excluded}  Total: {kept+excluded}", file=sys.stderr)
