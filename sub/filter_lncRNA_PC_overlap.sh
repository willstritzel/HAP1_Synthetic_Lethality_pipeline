#!/bin/bash
# filter_lncRNA_PC_overlap.sh
#
# Filters lncRNA_introns_final.bed to remove all entries for lncRNA genes
# that overlap protein-coding gene bodies (strand-ignored: removes both
# sense and antisense overlaps).
#
# Prerequisites:
#   - gencode.v49.basic.annotation.gtf.gz must already be downloaded to ANNOT_DIR
#   - gencode.v49.long_noncoding_RNAs.gtf must be present in ANNOT_DIR
#   - bedtools available in hap1 conda env
#   - python3 available in hap1 conda env
#
# Run from any directory; all outputs go to ANNOT_DIR.

set -euo pipefail

ANNOT_DIR="/scratch/alpine/wist9668/reference/annotation/lnc_RNA"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PC_GTF="${ANNOT_DIR}/gencode.v49.basic.annotation.gtf.gz"
LNCRNA_GTF="${ANNOT_DIR}/gencode.v49.long_noncoding_RNAs.gtf"
INPUT_BED="${ANNOT_DIR}/lncRNA_introns_final.bed"

PC_BED="${ANNOT_DIR}/protein_coding_genes_v49_basic.bed"
LNCRNA_GENES_BED="${ANNOT_DIR}/lncRNA_genes_v49.bed"
OVERLAP_LIST="${ANNOT_DIR}/lncRNA_genes_overlapping_PC.txt"
OUTPUT_BED="${ANNOT_DIR}/lncRNA_introns_no_PC_overlap.bed"

ml miniforge

# --- Step 1: Extract protein-coding gene bodies from GENCODE v49 basic GTF ---
echo "[1/4] Extracting protein-coding gene bodies..."
zcat "${PC_GTF}" \
  | awk '$3=="gene" && /gene_type "protein_coding"/' \
  | awk 'BEGIN{OFS="\t"} {
      match($0, /gene_name "([^"]+)"/, a);
      print $1, $4-1, $5, a[1], ".", $7
    }' \
  | sort -k1,1 -k2,2n \
  > "${PC_BED}"
echo "  Protein-coding genes: $(wc -l < "${PC_BED}")"

# --- Step 2: Extract lncRNA gene bodies from GENCODE v49 lncRNA GTF ---
echo "[2/4] Extracting lncRNA gene bodies..."
awk '$3=="gene"' "${LNCRNA_GTF}" \
  | awk 'BEGIN{OFS="\t"} {
      match($0, /gene_name "([^"]+)"/, a);
      print $1, $4-1, $5, a[1], ".", $7
    }' \
  | sort -k1,1 -k2,2n \
  > "${LNCRNA_GENES_BED}"
echo "  lncRNA genes: $(wc -l < "${LNCRNA_GENES_BED}")"

# --- Step 3: Find lncRNA genes overlapping any protein-coding gene body ---
echo "[3/4] Finding overlapping lncRNA genes (strand-ignored)..."
conda run -n hap1 bedtools intersect \
  -a "${LNCRNA_GENES_BED}" \
  -b "${PC_BED}" \
  -u \
  | cut -f4 | sort -u > "${OVERLAP_LIST}"
echo "  lncRNA genes overlapping PC genes: $(wc -l < "${OVERLAP_LIST}")"

# --- Step 4: Filter lncRNA introns BED ---
echo "[4/4] Filtering lncRNA_introns_final.bed..."
conda run -n hap1 python3 "${SCRIPT_DIR}/filter_lncRNA_introns.py" \
  "${OVERLAP_LIST}" \
  "${INPUT_BED}" \
  "${OUTPUT_BED}"

echo ""
echo "Done. Output: ${OUTPUT_BED}"
echo "  Input entries:  $(wc -l < "${INPUT_BED}")"
echo "  Output entries: $(wc -l < "${OUTPUT_BED}")"
