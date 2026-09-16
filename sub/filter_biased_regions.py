#!/usr/bin/env python3
"""
filter_biased_regions.py

Converts strand-biased regions from CM chromosome coordinates to UCSC chr names,
annotates against genomic features using bedtools, and filters based on user criteria.

Default behavior: protein-coding gene overlap is annotated and regions overlapping
PC genes are excluded. lncRNA genes are annotated as an informational column only.

Additional annotations can be added with --annotation LABEL:PATH.
Filtering logic can be relaxed (--no-exclude-pc) or extended (--exclude-lnc,
--exclude-all-annotated, --min-distance).
"""

import argparse
import os
import subprocess
import sys
import tempfile

DEFAULT_MAPPING = (
    "/scratch/alpine/wist9668/reference/annotation/lnc_RNA/gencode_to_cm_mapping.txt"
)
DEFAULT_PC_BED = (
    "/scratch/alpine/wist9668/reference/annotation/lnc_RNA/protein_coding_genes_v49_basic.bed"
)
DEFAULT_PC_GTF = (
    "/scratch/alpine/wist9668/reference/annotation/lnc_RNA/gencode.v49.basic.annotation.gtf.gz"
)
DEFAULT_LNC_BED = (
    "/scratch/alpine/wist9668/reference/annotation/lnc_RNA/lncRNA_genes_v49.bed"
)


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, help="BED from detect_strand_bias.py")
    p.add_argument("--output", required=True, help="Output filtered annotated BED")
    p.add_argument("--summary", required=True, help="TSV summary report path")
    p.add_argument(
        "--mapping",
        default=DEFAULT_MAPPING,
        help="gencode_to_cm_mapping.txt (chr<TAB>CM format) (default: %(default)s)",
    )
    p.add_argument(
        "--annotation",
        action="append",
        metavar="LABEL:PATH",
        type=lambda s: s.split(":", 1),
        help=(
            "Additional annotation BED in chr space with label (repeatable). "
            "Example: --annotation enh:/path/to/enhancers.bed"
        ),
    )
    p.add_argument(
        "--no-exclude-pc",
        action="store_true",
        default=False,
        help="Disable default exclusion of protein-coding gene overlaps",
    )
    p.add_argument(
        "--exclude-lnc",
        action="store_true",
        default=False,
        help="Exclude regions overlapping lncRNA gene bodies",
    )
    p.add_argument(
        "--exclude-all-annotated",
        action="store_true",
        default=False,
        help="Exclude regions overlapping any provided annotation",
    )
    p.add_argument(
        "--min-distance",
        type=int,
        default=0,
        help=(
            "Require regions to be >= N bp from any annotation "
            "(default: 0 = no distance filter)"
        ),
    )
    p.add_argument(
        "--keep-unplaced",
        action="store_true",
        default=False,
        help="Include alt contig regions in output (default: drop them)",
    )
    p.add_argument(
        "--bedtools",
        default="bedtools",
        help="Path to bedtools binary (default: bedtools on PATH)",
    )
    return p.parse_args()


def load_chr_mapping(mapping_path):
    """
    Read gencode_to_cm_mapping.txt (format: chr<TAB>CM).
    Returns (chr_to_cm, cm_to_chr).
    """
    chr_to_cm = {}
    cm_to_chr = {}
    with open(mapping_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            chr_name, cm_name = parts[0], parts[1]
            chr_to_cm[chr_name] = cm_name
            cm_to_chr[cm_name] = chr_name
    return chr_to_cm, cm_to_chr


def read_biased_regions(bed_path):
    """
    Read the 9-column output of detect_strand_bias.py.
    Returns list of dicts with keys: chrom, start, end, direction, sense_count,
    antisense_count, sense_fraction, p_value, fdr. Skips '#' header lines.
    """
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            regions.append({
                "chrom": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "direction": parts[3],
                "sense_count": int(parts[4]),
                "antisense_count": int(parts[5]),
                "sense_fraction": float(parts[6]),
                "p_value": float(parts[7]),
                "fdr": float(parts[8]),
            })
    return regions


def convert_to_chr(regions, cm_to_chr, keep_unplaced):
    """
    Replace CM accession names with UCSC chr names.
    Returns (converted_regions, list_of_dropped_cm_names).
    """
    converted = []
    dropped_cm = []
    for r in regions:
        chr_name = cm_to_chr.get(r["chrom"])
        if chr_name is None:
            dropped_cm.append(r["chrom"])
            if keep_unplaced:
                converted.append(dict(r))
        else:
            converted.append({**r, "chrom": chr_name})
    return converted, dropped_cm


def write_temp_bed9(regions, suffix=".bed"):
    """Write the 9-column biased regions to a temp file. Returns path."""
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=suffix, delete=False, prefix="sbias_"
    )
    for r in regions:
        tmp.write(
            f"{r['chrom']}\t{r['start']}\t{r['end']}\t"
            f"{r['direction']}\t{r['sense_count']}\t{r['antisense_count']}\t"
            f"{r['sense_fraction']}\t{r['p_value']}\t{r['fdr']}\n"
        )
    tmp.close()
    return tmp.name


def sort_bed(bed_path, bedtools):
    """Sort a BED file with bedtools sort. Returns path to a new sorted temp file."""
    sorted_tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".sorted.bed", delete=False, prefix="sbias_sorted_"
    )
    sorted_tmp.close()
    result = subprocess.run(
        [bedtools, "sort", "-i", bed_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bedtools sort failed:\n{result.stderr}")
    with open(sorted_tmp.name, "w") as out:
        out.write(result.stdout)
    return sorted_tmp.name


def get_overlapping_keys(query_bed, subject_bed, bedtools):
    """
    Run bedtools intersect -u to find query entries overlapping subject.
    Returns set of (chrom, start, end) tuples.
    """
    result = subprocess.run(
        [bedtools, "intersect", "-a", query_bed, "-b", subject_bed, "-u"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bedtools intersect failed:\n{result.stderr}")
    overlapping = set()
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        overlapping.add((parts[0], int(parts[1]), int(parts[2])))
    return overlapping


def get_closest_distances(query_bed, subject_bed, bedtools):
    """
    Run bedtools closest -d to get the distance from each query entry to the
    nearest subject feature. Returns dict (chrom, start, end) -> int distance.
    Distance = 0 means overlap; -1 means no features on the chromosome.
    """
    result = subprocess.run(
        [bedtools, "closest", "-a", query_bed, "-b", subject_bed, "-d"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bedtools closest failed:\n{result.stderr}")
    distances = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        key = (parts[0], int(parts[1]), int(parts[2]))
        distances[key] = int(parts[-1])
    return distances


def ensure_pc_bed(pc_bed_path, pc_gtf_path):
    """
    Return pc_bed_path if it exists and is non-empty.
    Otherwise extract protein-coding gene bodies from the GENCODE GTF on-the-fly.
    """
    if os.path.isfile(pc_bed_path) and os.path.getsize(pc_bed_path) > 0:
        return pc_bed_path

    print(
        f"  PC gene BED not found; extracting from {pc_gtf_path} ...",
        file=sys.stderr,
    )
    cmd = (
        f'zcat "{pc_gtf_path}" '
        r'| awk \'$3=="gene" && /gene_type "protein_coding"/\' '
        r'| awk \'BEGIN{OFS="\t"} {match($0, /gene_name "([^"]+)"/, a); print $1, $4-1, $5, a[1], ".", $7}\' '
        f'| sort -k1,1 -k2,2n > "{pc_bed_path}"'
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"PC gene extraction failed:\n{result.stderr}")
    return pc_bed_path


def annotate_regions(regions, annotations, bedtools, min_distance):
    """
    For each (label, bed_path) in annotations, run bedtools to determine:
      - r['annot_{label}']: bool, whether the region overlaps the annotation
      - r['dist_{label}']:  int, distance to nearest feature (only if min_distance > 0)
    Mutates region dicts in place.
    """
    if not regions:
        return

    query_tmp = write_temp_bed9(regions)
    sorted_query = sort_bed(query_tmp, bedtools)
    os.unlink(query_tmp)

    try:
        for label, annot_bed in annotations:
            print(f"  [{label}] {annot_bed}", file=sys.stderr)
            overlapping = get_overlapping_keys(sorted_query, annot_bed, bedtools)
            distances = (
                get_closest_distances(sorted_query, annot_bed, bedtools)
                if min_distance > 0
                else {}
            )
            for r in regions:
                key = (r["chrom"], r["start"], r["end"])
                r[f"annot_{label}"] = key in overlapping
                r[f"dist_{label}"] = distances.get(key) if min_distance > 0 else None
    finally:
        os.unlink(sorted_query)


def apply_filters(regions, exclude_pc, exclude_lnc, exclude_all_annotated,
                  min_distance, annotation_labels):
    """Return the subset of regions that pass all active filters."""
    filtered = []
    for r in regions:
        keep = True

        if exclude_all_annotated:
            if any(r.get(f"annot_{lbl}", False) for lbl in annotation_labels):
                keep = False
        else:
            if exclude_pc and r.get("annot_pc", False):
                keep = False
            if exclude_lnc and r.get("annot_lnc", False):
                keep = False

        if keep and min_distance > 0:
            for lbl in annotation_labels:
                dist = r.get(f"dist_{lbl}")
                if dist is not None and 0 <= dist < min_distance:
                    keep = False
                    break

        if keep:
            filtered.append(r)
    return filtered


def write_filtered_bed(regions, outpath, annotation_labels):
    """Write annotated, filtered BED with comment header."""
    base = [
        "chrom", "start", "end", "direction", "sense_count",
        "antisense_count", "sense_fraction", "p_value", "fdr",
    ]
    annot_cols = []
    for lbl in annotation_labels:
        annot_cols += [f"annot_{lbl}", f"dist_{lbl}"]

    with open(outpath, "w") as out:
        out.write("#" + "\t".join(base + annot_cols) + "\n")
        for r in regions:
            row = [str(r.get(c, "")) for c in base]
            for lbl in annotation_labels:
                row.append("True" if r.get(f"annot_{lbl}", False) else "False")
                dist = r.get(f"dist_{lbl}")
                row.append(str(dist) if dist is not None else "")
            out.write("\t".join(row) + "\n")


def write_summary(all_regions, filtered_regions, annotation_labels, outpath):
    """Write a TSV summary of region counts and annotation overlaps."""
    with open(outpath, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_regions_input\t{len(all_regions)}\n")
        out.write(f"total_regions_output\t{len(filtered_regions)}\n")
        out.write(f"regions_excluded\t{len(all_regions) - len(filtered_regions)}\n")

        for direction in ("sense_bias", "antisense_bias"):
            n_in = sum(1 for r in all_regions if r["direction"] == direction)
            n_out = sum(1 for r in filtered_regions if r["direction"] == direction)
            out.write(f"input_{direction}\t{n_in}\n")
            out.write(f"output_{direction}\t{n_out}\n")

        for lbl in annotation_labels:
            n_ov = sum(1 for r in all_regions if r.get(f"annot_{lbl}", False))
            pct = f"{100 * n_ov / len(all_regions):.1f}" if all_regions else "NA"
            out.write(f"overlap_with_{lbl}\t{n_ov}\n")
            out.write(f"pct_overlap_with_{lbl}\t{pct}\n")


def main():
    args = parse_args()

    print(f"[filter_biased_regions] Input:   {args.input}", file=sys.stderr)
    print(f"[filter_biased_regions] Output:  {args.output}", file=sys.stderr)
    print(f"[filter_biased_regions] Summary: {args.summary}", file=sys.stderr)

    _, cm_to_chr = load_chr_mapping(args.mapping)

    regions = read_biased_regions(args.input)
    print(f"[filter_biased_regions] {len(regions):,} regions loaded", file=sys.stderr)

    regions, dropped = convert_to_chr(regions, cm_to_chr, args.keep_unplaced)
    if dropped:
        print(
            f"  {len(set(dropped))} unplaced contig(s) dropped "
            f"({len(dropped)} region(s))",
            file=sys.stderr,
        )
    print(
        f"  {len(regions):,} regions retained after coordinate conversion",
        file=sys.stderr,
    )

    # Build annotation list: PC always first, lncRNA second if available
    annotations = list(args.annotation) if args.annotation else []
    present_labels = [lbl for lbl, _ in annotations]

    if "pc" not in present_labels:
        pc_bed = ensure_pc_bed(DEFAULT_PC_BED, DEFAULT_PC_GTF)
        annotations.insert(0, ["pc", pc_bed])

    if "lnc" not in present_labels and os.path.isfile(DEFAULT_LNC_BED):
        annotations.append(["lnc", DEFAULT_LNC_BED])

    annotation_labels = [lbl for lbl, _ in annotations]
    exclude_pc = not args.no_exclude_pc

    print(
        f"[filter_biased_regions] Annotating with bedtools "
        f"(labels: {', '.join(annotation_labels)})...",
        file=sys.stderr,
    )
    annotate_regions(regions, annotations, args.bedtools, args.min_distance)

    filtered = apply_filters(
        regions,
        exclude_pc=exclude_pc,
        exclude_lnc=args.exclude_lnc,
        exclude_all_annotated=args.exclude_all_annotated,
        min_distance=args.min_distance,
        annotation_labels=annotation_labels,
    )
    print(
        f"[filter_biased_regions] {len(filtered):,} / {len(regions):,} regions pass filters",
        file=sys.stderr,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    write_filtered_bed(filtered, args.output, annotation_labels)
    write_summary(regions, filtered, annotation_labels, args.summary)

    print(f"[filter_biased_regions] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
