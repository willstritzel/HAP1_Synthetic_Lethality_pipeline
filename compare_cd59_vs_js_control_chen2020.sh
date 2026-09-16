#!/bin/bash
# Comparative haploid screen analysis: CD59 (experimental) vs js_hap1_control (unselected)
# Using Chen et al. 2020 ORF annotation (introns + exons)
# Following Davis et al. 2015 (Cell Reports doi:10.1016/j.celrep.2015.05.026)
#
# Run batch_job_cd59_chen2020.sh and batch_job_ctrl_chen2020.sh first and wait
# for both to complete before running this script.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/sub"

REF_DIR="/scratch/alpine/wist9668/reference/annotation/chen2020"
EXONS_BED="${REF_DIR}/chen2020_exon_CM-mapped.bed"

DATA_DIR="/scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens"
CD59_DIR="${DATA_DIR}/JS_HAP1_CD59_chen2020_output/replicate_1"
CTRL_DIR="${DATA_DIR}/JS_ControlData_chen2020_output/replicate_1"

CD59_BWT="${CD59_DIR}/replicate.mapped_unique.bwt"
CTRL_BWT="${CTRL_DIR}/replicate.mapped_unique.bwt"

CD59_INTRON="${CD59_DIR}/chen2020/sense_vs_antisense_counts.csv"
CTRL_INTRON="${CTRL_DIR}/chen2020/sense_vs_antisense_counts.csv"

OUT_DIR="${DATA_DIR}/CD59_vs_JScontrol_chen2020_comparative_screen"
mkdir -p "${OUT_DIR}"

CD59_EXON="${OUT_DIR}/CD59_exon_counts.tsv"
CTRL_EXON="${OUT_DIR}/control_exon_counts.tsv"
RESULTS="${OUT_DIR}/CD59_vs_JScontrol_chen2020_results.tsv"

for f in "${EXONS_BED}" "${CD59_BWT}" "${CTRL_BWT}" "${CD59_INTRON}" "${CTRL_INTRON}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required file not found: ${f}"
        exit 1
    fi
done

echo "[1/3] Counting exonic insertions for CD59 screen ..."
ml miniforge && conda run -n hap1 python3 "${SCRIPT_DIR}/count_exon_insertions.py" \
    --bwt   "${CD59_BWT}" \
    --exons "${EXONS_BED}" \
    --out   "${CD59_EXON}"

echo "[1/3] Counting exonic insertions for control ..."
ml miniforge && conda run -n hap1 python3 "${SCRIPT_DIR}/count_exon_insertions.py" \
    --bwt   "${CTRL_BWT}" \
    --exons "${EXONS_BED}" \
    --out   "${CTRL_EXON}"

echo "[2/3] Running comparative screen analysis ..."
ml miniforge && conda run -n hap1 python3 "${SCRIPT_DIR}/compare_experimental_vs_control.py" \
    --intron-exp  "${CD59_INTRON}" \
    --intron-ctrl "${CTRL_INTRON}" \
    --exon-exp    "${CD59_EXON}" \
    --exon-ctrl   "${CTRL_EXON}" \
    --out         "${RESULTS}" \
    --pthresh     1e-5

echo ""
echo "Results: ${RESULTS}"
