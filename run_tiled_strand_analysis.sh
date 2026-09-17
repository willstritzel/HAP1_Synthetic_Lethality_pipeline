#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=tiled_strand_bias
#SBATCH --output=logs/tiled_strand_bias.%j.out
#SBATCH --time=006:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

# run_tiled_strand_analysis.sh
#
# Aggregates raw insertion alignments from tiled HAP1 screen replicates,
# detects strand-biased genomic regions via sliding window binomial test,
# and filters out canonical genomic annotations.
#
# Usage:
#   sbatch run_tiled_strand_analysis.sh [options]
#   bash   run_tiled_strand_analysis.sh [options]   # local test run
#
# Options (key=value format):
#   --screen-name=STR      Label for output files (default: hap1-4_tiled)
#   --rep-dir=PATH         Directory with replicate_{1..4}/ subdirs
#   --outdir=PATH          Output directory (default: see FINAL_DIR below)
#   --window-size=INT      Sliding window size in bp (default: 10000)
#   --step-size=INT        Step size in bp (default: 1000)
#   --min-insertions=INT   Min insertions per window to test (default: 10)
#   --fdr-threshold=FLOAT  BH FDR threshold (default: 0.05)
#   --merge-gap=INT        Max gap to merge adjacent significant windows (default: 1000)
#   --no-exclude-pc        Do not filter out protein-coding gene overlaps
#   --exclude-lnc          Also filter out lncRNA gene overlaps
#   --exclude-all          Filter out all annotated regions
#   --min-distance=INT     Require N bp distance from any annotation (default: 0)

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
SCREEN_NAME="hap1-4_tiled"
REP_DIR="/scratch/alpine/wist9668/haploid_scratch/outputs/tiled_2kb_hap1-4_output"
FINAL_DIR="/scratch/alpine/wist9668/haploid_scratch/outputs"
WINDOW_SIZE=10000
STEP_SIZE=1000
MIN_INSERTIONS=10
FDR_THRESHOLD=0.05
MERGE_GAP=1000
FILTER_FLAGS=""
# ─────────────────────────────────────────────────────────────────────────────

# ── Parse arguments ───────────────────────────────────────────────────────────
for i in "$@"; do
    case $i in
        --screen-name=*)    SCREEN_NAME="${i#*=}" ;;
        --rep-dir=*)        REP_DIR="${i#*=}" ;;
        --outdir=*)         FINAL_DIR="${i#*=}" ;;
        --window-size=*)    WINDOW_SIZE="${i#*=}" ;;
        --step-size=*)      STEP_SIZE="${i#*=}" ;;
        --min-insertions=*) MIN_INSERTIONS="${i#*=}" ;;
        --fdr-threshold=*)  FDR_THRESHOLD="${i#*=}" ;;
        --merge-gap=*)      MERGE_GAP="${i#*=}" ;;
        --no-exclude-pc)    FILTER_FLAGS="${FILTER_FLAGS} --no-exclude-pc" ;;
        --exclude-lnc)      FILTER_FLAGS="${FILTER_FLAGS} --exclude-lnc" ;;
        --exclude-all)      FILTER_FLAGS="${FILTER_FLAGS} --exclude-all-annotated" ;;
        --min-distance=*)   FILTER_FLAGS="${FILTER_FLAGS} --min-distance=${i#*=}" ;;
        *) echo "Unknown option: $i" >&2; exit 1 ;;
    esac
done

# ── Setup ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
OUTDIR="${FINAL_DIR}/${SCREEN_NAME}_strand_bias"
mkdir -p "${OUTDIR}" logs

INSERTIONS_BED="${OUTDIR}/${SCREEN_NAME}_insertions.bed"
BIAS_BED="${OUTDIR}/${SCREEN_NAME}_strand_bias_raw.bed"
FILTERED_BED="${OUTDIR}/${SCREEN_NAME}_biased_regions.bed"
SUMMARY_TSV="${OUTDIR}/${SCREEN_NAME}_biased_regions_summary.tsv"

echo "========================================================"
echo " Tiled strand bias analysis"
echo " Screen:      ${SCREEN_NAME}"
echo " Rep dir:     ${REP_DIR}"
echo " Output dir:  ${OUTDIR}"
echo " Window:      ${WINDOW_SIZE} bp / ${STEP_SIZE} bp step"
echo " Min inserts: ${MIN_INSERTIONS}"
echo " FDR:         ${FDR_THRESHOLD}"
echo " Merge gap:   ${MERGE_GAP} bp"
echo " Started:     $(date)"
echo "========================================================"

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh

# ── Step 1: Aggregate and de-duplicate insertions ─────────────────────────────
echo ""
echo "[$(date '+%H:%M:%S')] Step 1: Aggregating insertions from replicates..."
conda run -n hap1 python3 "${SCRIPT_DIR}/sub/aggregate_insertions.py" \
    --rep-dir="${REP_DIR}" \
    --output="${INSERTIONS_BED}" \
    --screen-name="${SCREEN_NAME}"

echo "[$(date '+%H:%M:%S')] Step 1 complete. Lines: $(wc -l < "${INSERTIONS_BED}")"

# ── Step 2: Detect strand-biased regions ──────────────────────────────────────
echo ""
echo "[$(date '+%H:%M:%S')] Step 2: Detecting strand-biased regions..."
conda run -n hap1 python3 "${SCRIPT_DIR}/sub/detect_strand_bias.py" \
    --input="${INSERTIONS_BED}" \
    --output="${BIAS_BED}" \
    --window-size="${WINDOW_SIZE}" \
    --step-size="${STEP_SIZE}" \
    --min-insertions="${MIN_INSERTIONS}" \
    --fdr-threshold="${FDR_THRESHOLD}" \
    --merge-gap="${MERGE_GAP}"

BIAS_COUNT=$(grep -v "^#" "${BIAS_BED}" | wc -l)
echo "[$(date '+%H:%M:%S')] Step 2 complete. Biased regions (pre-filter): ${BIAS_COUNT}"

# ── Step 3: Annotate and filter ───────────────────────────────────────────────
echo ""
echo "[$(date '+%H:%M:%S')] Step 3: Annotating and filtering biased regions..."
conda run -n hap1 python3 "${SCRIPT_DIR}/sub/filter_biased_regions.py" \
    --input="${BIAS_BED}" \
    --output="${FILTERED_BED}" \
    --summary="${SUMMARY_TSV}" \
    ${FILTER_FLAGS}

FINAL_COUNT=$(grep -v "^#" "${FILTERED_BED}" | wc -l)
echo "[$(date '+%H:%M:%S')] Step 3 complete. Final regions: ${FINAL_COUNT}"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "========================================================"
echo " Analysis complete: $(date)"
echo " Output directory:  ${OUTDIR}"
echo "   Insertions BED:  $(basename "${INSERTIONS_BED}")"
echo "   Raw bias BED:    $(basename "${BIAS_BED}")"
echo "   Filtered BED:    $(basename "${FILTERED_BED}")"
echo "   Summary TSV:     $(basename "${SUMMARY_TSV}")"
echo " Regions (pre-filter):  ${BIAS_COUNT}"
echo " Regions (final):       ${FINAL_COUNT}"
echo "========================================================"
