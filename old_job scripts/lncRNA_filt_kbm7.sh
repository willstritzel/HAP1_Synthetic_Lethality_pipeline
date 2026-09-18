#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=lncRNA_filt_kbm7
#SBATCH --output=logs/lncRNA_filt_kbm7.%j.out
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# ── Job-specific settings ────────────────────────────────────────────────────
SCREEN_NAME=lncRNA_filt_kbm7
GENE_REF=/scratch/alpine/wist9668/haploid_scratch/reference/annotation/lnc_RNA/lncRNA_introns_no_PC_overlap_CM.bed
REF_ID=lncRNA_noPC
# ─────────────────────────────────────────────────────────────────────────────

# Clean up previous runs
rm -rf /scratch/alpine/wist9668/haploid_scratch/tmp/${SCREEN_NAME}_output/
rm -rf /scratch/alpine/wist9668/haploid_scratch/outputs/${SCREEN_NAME}_output/

# Array of input files mapped to replicate numbers
FILES=(
    "/scratch/alpine/wist9668/KBM7/SRR2047137.fastq.gz"
    "/scratch/alpine/wist9668/KBM7/SRR2047138.fastq.gz"
    "/scratch/alpine/wist9668/KBM7/SRR2047149.fastq.gz"
)

# Run each replicate
for REP in {1..3}; do
    echo "==============================================="
    echo "Starting replicate ${REP} at $(date)"
    echo "==============================================="

    ./analyze_sli.sh -C -R=${REP} -S=${FILES[$((REP-1))]} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}

    echo "Completed replicate ${REP} at $(date)"
    echo ""
done

echo "All replicates completed at $(date)"
