#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=JOBNAME
#SBATCH --output=logs/JOBNAME.%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

#####################################################################
# Template for a 4-replicate screen. Copy, then edit the three
# job-specific settings and the FILES array.
#
# Cluster facts this header encodes (verified 2026-09-16):
#   partition acpu accepts only qos cpu-normal or cpu-long; the old
#   amilan/normal pair no longer exists. Short tests go to
#   --partition=atesting --qos=testing (16 CPUs, 1 h, 1 job).
#   There is no miniforge module and no mamba on PATH; conda must be
#   sourced from /curc/sw/anaconda3/2023.09.
#####################################################################

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# ── Job-specific settings ────────────────────────────────────────────────────
SCREEN_NAME=SCREENNAME
GENE_REF=/pl/active/ShenLab_PL/haploid_genetics/reference/annotation/ANNOTATION.bed
REF_ID=ref_id_label
# ─────────────────────────────────────────────────────────────────────────────

TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)

# Inputs: keep these on the PetaLibrary. Scratch copies get purged.
FILES=(
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/DATASET/SRR_REP1.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/DATASET/SRR_REP2.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/DATASET/SRR_REP3.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/DATASET/SRR_REP4.fastq.gz"
)

# ── Preflight ────────────────────────────────────────────────────────────────
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
