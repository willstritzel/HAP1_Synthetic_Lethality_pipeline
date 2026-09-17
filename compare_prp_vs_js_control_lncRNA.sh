#!/bin/bash
# Comparative haploid screen analysis: PrP (experimental) vs js_hap1_control (unselected)
# Using Chothani et al. 2025 lncRNA annotation (introns + exons)
# Following Davis et al. 2015 (Cell Reports doi:10.1016/j.celrep.2015.05.026)
#
# Run batch_job_prp_lncrna.sh and batch_job_ctrl_lncrna.sh first and wait
# for both to complete before running this script.

set -euo pipefail

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/sub"

REF_DIR="/scratch/alpine/wist9668/haploid_scratch/reference/annotation/lnc_RNA"
EXONS_BED="${REF_DIR}/chothani_lncRNA_exon_CM-mapped.bed"
INTRONS_BED="${REF_DIR}/chothani_lncRNA_intron_real_CM-mapped.bed"

DATA_DIR="/scratch/alpine/wist9668/haploid_scratch/screens"
PrP_DIR="${DATA_DIR}/JS_HAP1_PrP_lncRNA_output/replicate_1"
CTRL_DIR="${DATA_DIR}/JS_ControlData_lncRNA_output/replicate_1"

PrP_BWT="${PrP_DIR}/replicate.mapped_unique.bwt"
CTRL_BWT="${CTRL_DIR}/replicate.mapped_unique.bwt"

PrP_INTRON="${PrP_DIR}/chothani_lncRNA/sense_vs_antisense_counts.csv"
CTRL_INTRON="${CTRL_DIR}/chothani_lncRNA/sense_vs_antisense_counts.csv"

OUT_DIR="${DATA_DIR}/PrP_vs_JScontrol_lncRNA_comparative_screen"
mkdir -p "${OUT_DIR}"

PrP_EXON="${OUT_DIR}/PrP_exon_counts.tsv"
CTRL_EXON="${OUT_DIR}/control_exon_counts.tsv"
RESULTS="${OUT_DIR}/PrP_vs_JScontrol_lncRNA_results.tsv"

for f in "${EXONS_BED}" "${PrP_BWT}" "${CTRL_BWT}" "${PrP_INTRON}" "${CTRL_INTRON}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required file not found: ${f}"
        exit 1
    fi
done

echo "[1/3] Counting exonic insertions for PrP screen ..."
conda run -n hap1 python3 "${SCRIPT_DIR}/count_exon_insertions.py" \
    --bwt   "${PrP_BWT}" \
    --exons "${EXONS_BED}" \
    --out   "${PrP_EXON}"

echo "[1/3] Counting exonic insertions for control ..."
conda run -n hap1 python3 "${SCRIPT_DIR}/count_exon_insertions.py" \
    --bwt   "${CTRL_BWT}" \
    --exons "${EXONS_BED}" \
    --out   "${CTRL_EXON}"

echo "[2/3] Running comparative screen analysis ..."
conda run -n hap1 python3 "${SCRIPT_DIR}/compare_experimental_vs_control.py" \
    --intron-exp  "${PrP_INTRON}" \
    --intron-ctrl "${CTRL_INTRON}" \
    --exon-exp    "${PrP_EXON}" \
    --exon-ctrl   "${CTRL_EXON}" \
    --out         "${RESULTS}" \
    --pthresh     1e-5

echo ""
echo "Results: ${RESULTS}"
