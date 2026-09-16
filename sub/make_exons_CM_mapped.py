#!/usr/bin/env python3
"""
Convert exons.bed from NC chromosome naming to CM accession naming.

Reads the GRCh38 assembly report to build the NC→CM mapping, then remaps
column 1 of exons.bed. Rows on contigs not present in the mapping are dropped.

Usage:
    python3 make_exons_CM_mapped.py <assembly_report> <exons_bed> <output_bed>
"""
import sys

if len(sys.argv) != 4:
    print("Usage: make_exons_CM_mapped.py <assembly_report> <exons_bed> <output_bed>",
          file=sys.stderr)
    sys.exit(1)

assembly_report, exons_bed, output_bed = sys.argv[1], sys.argv[2], sys.argv[3]

# Build NC → CM mapping from assembly report
# Non-comment lines: col 4 = CM accession, col 6 = NC accession
nc_to_cm = {}
with open(assembly_report) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        cm = parts[4]
        nc = parts[6]
        if nc.startswith('NC_') and (cm.startswith('CM') or cm.startswith('NC')):
            nc_to_cm[nc] = cm

kept = dropped = 0
with open(exons_bed) as fin, open(output_bed, 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        parts = line.strip().split('\t')
        if not parts:
            continue
        nc = parts[0]
        cm = nc_to_cm.get(nc)
        if cm is None:
            dropped += 1
            continue
        parts[0] = cm
        fout.write('\t'.join(parts) + '\n')
        kept += 1

print(f"Done: {kept} exon entries kept, {dropped} dropped (unmapped contigs)",
      file=sys.stderr)
