#!/bin/bash
# Comparative haploid screen analysis: PrP (experimental) vs js_hap1_control (unselected)
# Following Davis et al. 2015 (Cell Reports doi:10.1016/j.celrep.2015.05.026)
#
# Independent of the existing strand-bias pipeline — does not read or modify
# any existing pipeline scripts or config files.

set -euo pipefail

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh

# ── Paths ──────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/sub"
PYTHON="conda run -n hap1 python3"

REF_DIR="/scratch/alpine/wist9668/haploid_scratch/reference/annotation"
ASSEMBLY_REPORT="${REF_DIR}/GCA_000001405.15_GRCh38_assembly_report.txt"
EXONS_NC="${REF_DIR}/exons.bed"
EXONS_CM="${REF_DIR}/exons_CM_mapped.bed"

DATA_DIR="/scratch/alpine/wist9668/haploid_scratch/outputs"
PrP_DIR="${DATA_DIR}/JS_HAP1_PrP_output/replicate_1"
CTRL_DIR="${DATA_DIR}/JS_ControlData-HAP1_output/replicate_1"

PrP_BWT="${PrP_DIR}/replicate.mapped_unique.bwt"
CTRL_BWT="${CTRL_DIR}/replicate.mapped_unique.bwt"

PrP_INTRON="${PrP_DIR}/PrP/sense_vs_antisense_counts.csv"
CTRL_INTRON="${CTRL_DIR}/GRCh38-NCBI_RefSeq/sense_vs_antisense_counts.csv"

OUT_DIR="${DATA_DIR}/PrP_vs_JScontrol_comparative_screen"
mkdir -p "${OUT_DIR}"

PrP_EXON="${OUT_DIR}/PrP_exon_counts.tsv"
CTRL_EXON="${OUT_DIR}/control_exon_counts.tsv"
RESULTS="${OUT_DIR}/PrP_vs_JScontrol_results.tsv"

# ── Step 1: Create exons_CM_mapped.bed (if not already done) ──────────────────

if [[ ! -f "${EXONS_CM}" ]]; then
    echo "[1/3] Creating exons_CM_mapped.bed ..."
    conda run -n hap1 python3 "${SCRIPT_DIR}/make_exons_CM_mapped.py" \
        "${ASSEMBLY_REPORT}" "${EXONS_NC}" "${EXONS_CM}"
else
    echo "[1/3] exons_CM_mapped.bed already exists, skipping."
fi

# ── Step 2: Count exonic insertions per gene ───────────────────────────────────

echo "[2/3] Counting exonic insertions for PrP screen ..."
conda run -n hap1 python3 "${SCRIPT_DIR}/count_exon_insertions.py" \
    --bwt   "${PrP_BWT}" \
    --exons "${EXONS_CM}" \
    --out   "${PrP_EXON}"

echo "[2/3] Counting exonic insertions for control ..."
conda run -n hap1 python3 "${SCRIPT_DIR}/count_exon_insertions.py" \
    --bwt   "${CTRL_BWT}" \
    --exons "${EXONS_CM}" \
    --out   "${CTRL_EXON}"

# ── Step 3: Run comparative Fisher test ───────────────────────────────────────

echo "[3/3] Running comparative screen analysis ..."
conda run -n hap1 python3 "${SCRIPT_DIR}/compare_experimental_vs_control.py" \
    --intron-exp  "${PrP_INTRON}" \
    --intron-ctrl "${CTRL_INTRON}" \
    --exon-exp    "${PrP_EXON}" \
    --exon-ctrl   "${CTRL_EXON}" \
    --out         "${RESULTS}" \
    --pthresh     1e-5

echo ""
echo "Results: ${RESULTS}"
