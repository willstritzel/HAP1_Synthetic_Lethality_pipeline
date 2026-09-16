#!/usr/bin/env python3
"""
Convert chromosome names in a BED file using a mapping file.
Entries on contigs not present in the mapping are dropped.

Usage: rename_chromosomes.py <mapping_file> <input_bed> <output_bed>

Mapping file format (tab-separated, 3 columns):
  <number>  <chr_name>  <accession>
  e.g.: 1   chr1   CM000663.2
"""
import sys

if len(sys.argv) != 4:
    print("Usage: rename_chromosomes.py <mapping_file> <input_bed> <output_bed>",
          file=sys.stderr)
    sys.exit(1)

mapping_file, input_bed, output_bed = sys.argv[1], sys.argv[2], sys.argv[3]

mapping = {}
with open(mapping_file) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            mapping[parts[0]] = parts[1]  # chr_name -> accession

kept = dropped = 0
with open(input_bed) as fin, open(output_bed, "w") as fout:
    for line in fin:
        chrom = line.split("\t")[0]
        if chrom in mapping:
            fout.write(mapping[chrom] + "\t" + line.split("\t", 1)[1])
            kept += 1
        else:
            dropped += 1

print(f"Kept: {kept}  Dropped (unmapped contig): {dropped}", file=sys.stderr)
