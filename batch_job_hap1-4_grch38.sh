#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=grch38_hap1-4
#SBATCH --output=logs/grch38_hap1-4.%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

#####################################################################
# Standard HAP1 gene-trap screen: the four unselected HAP1 replicates
# (Blomen et al. 2015, SRR2047158-161) against the genome-wide GRCh38
# intron annotation.
#
# This is the reference run of the pipeline - the sense/antisense
# readout on real gene models, as opposed to the tiled-window null or
# the chen2020 / lncRNA feature sets. The screen name ControlData-HAP1
# matches the four HAP1_control_* entries in sub/analytical_settings.py
# and the legacy archive at outputs/legacy/ControlData-HAP1_output.
#
# -C (control-creation mode) is correct here and is not a workaround:
# hap1-4 IS the unselected control dataset, so there is no control to
# compare against and normalize.py / test_replicate_vs_control.py are
# rightly skipped. See the F07 note in the project docs.
#####################################################################

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# ── Job-specific settings ────────────────────────────────────────────────────
SCREEN_NAME=ControlData-HAP1
GENE_REF=/pl/active/ShenLab_PL/haploid_genetics/reference/annotation/GRCh38_introns_CM_mapped.bed
REF_ID=GRCh38-NCBI_RefSeq
# ─────────────────────────────────────────────────────────────────────────────

# Output locations come from the settings file, so a change there does not
# have to be mirrored into every wrapper by hand.
TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)

# Array of input files mapped to replicate numbers
FILES=(
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047158.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047159.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047160.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047161.fastq.gz"
)

# ── Preflight ────────────────────────────────────────────────────────────────
# Fail before the first alignment rather than four hours in. The 2026 scratch
# purge left every wrapper pointing at annotation BEDs and FASTQs that no
# longer existed (F04, F05).
echo "Preflight checks:"
MISSING=0
if [[ -s $GENE_REF ]]; then
    echo "  gene_ref OK   ($(wc -l < "$GENE_REF") features) $GENE_REF"
else
    echo "  gene_ref MISSING OR EMPTY: $GENE_REF"; MISSING=1
fi
for f in "${FILES[@]}"; do
    if [[ -s $f ]]; then echo "  fastq OK      $f"
    else echo "  fastq MISSING: $f"; MISSING=1; fi
done
if [[ $MISSING -ne 0 ]]; then
    echo "Preflight failed - nothing submitted to bowtie. Exiting."; exit 1
fi
echo "  tmp_dir       $TMP_DIR"
echo "  final_dir     $FINAL_DIR"
mkdir -p "$TMP_DIR" "$FINAL_DIR"

# Clean up previous runs
rm -rf "${TMP_DIR}${SCREEN_NAME}_output/"
rm -rf "${FINAL_DIR}${SCREEN_NAME}_output/"

# Run each replicate
for REP in {1..4}; do
    echo "==============================================="
    echo "Starting replicate ${REP} at $(date)"
    echo "==============================================="

    ./analyze_sli.sh -C -R=${REP} -S=${FILES[$((REP-1))]} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}

    echo "Completed replicate ${REP} at $(date)"
    echo ""
done

echo "All replicates completed at $(date)"

# ── Archive to PetaLibrary ───────────────────────────────────────────────────
# final_dir is on scratch and will be purged. Do not remove this step.
bash sub/archive_screen.sh "${SCREEN_NAME}" "${REF_ID}" "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out"
