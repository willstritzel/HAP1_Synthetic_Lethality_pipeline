#!/usr/bin/env python3
"""
Extract exon intervals for lncRNA ORFs from the Chothani et al. 2025 catalogue
and write a CM-chromosome-named BED6 file, complementary to
chothani_lncRNA_intron_real_CM-mapped.bed.

Filters (matching the intron BED):
    orf_type == "lncRNA"  AND  orf_length (bp) > 33

Each row in the catalogue can have multiple exons (comma-separated starts/ends).
All exons for each passing ORF are output — one BED row per exon, strand-agnostic
deduplication is NOT applied here (bedtools intersect handles counting later).

Chromosome mapping:
    catalogue chrm column uses plain numeric names (1, 2, ... X, Y, MT)
    → mapped to CM accessions (CM000663.2, ...) via GRCh38 assembly report
    Rows whose chrm name has no CM accession (alt contigs etc.) are dropped.

Usage:
    python3 make_lncrna_exons_CM_mapped.py \\
        <assembly_report.txt> \\
        <"Chothani et al. 2025 catalogue.xlsx"> \\
        <output_exon_CM-mapped.bed>
"""
import sys
import pandas as pd


def load_chr_map(assembly_report):
    """Build {sequence_name: CM_accession} for assembled chromosomes only."""
    mapping = {}
    with open(assembly_report) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            if parts[1] != 'assembled-molecule':
                continue
            seq_name = parts[0]   # "1", "2", ..., "X", "Y", "MT"
            cm_acc   = parts[4]   # "CM000663.2", ...
            mapping[seq_name] = cm_acc
    return mapping


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    assembly_report, xlsx_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

    chr_map = load_chr_map(assembly_report)
    print(f"Chromosome map: {len(chr_map)} assembled sequences", flush=True)

    df = pd.read_excel(xlsx_path, sheet_name='Supplementary Table 2')
    print(f"Loaded {len(df):,} rows from catalogue", flush=True)

    filt = df[(df['orf_type'] == 'lncRNA') & (df['orf_length (bp)'] > 33)].copy()
    print(f"After filter (lncRNA, orf_length > 33): {len(filt):,} ORFs", flush=True)

    n_written = 0
    n_dropped_chr = 0
    n_dropped_parse = 0

    with open(out_path, 'w') as out:
        for _, row in filt.iterrows():
            chrm   = str(row['chrm']).strip()
            cm     = chr_map.get(chrm)
            if cm is None:
                n_dropped_chr += 1
                continue

            name   = str(row['releasev45_id']).strip()
            strand = str(row['strand']).strip()

            try:
                starts = [int(x) for x in str(row['starts (0-based)']).split(',')]
                ends   = [int(x) for x in str(row['ends (0-based)']).split(',')]
            except (ValueError, AttributeError):
                n_dropped_parse += 1
                continue

            if len(starts) != len(ends):
                n_dropped_parse += 1
                continue

            for s, e in zip(starts, ends):
                out.write(f"{cm}\t{s}\t{e}\t{name}\t.\t{strand}\n")
                n_written += 1

    print(f"Exon BED rows written : {n_written:,}", flush=True)
    print(f"ORFs dropped (no CM)  : {n_dropped_chr}", flush=True)
    print(f"ORFs dropped (parse)  : {n_dropped_parse}", flush=True)
    print(f"Output: {out_path}", flush=True)


if __name__ == '__main__':
    main()
