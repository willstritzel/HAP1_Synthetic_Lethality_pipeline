#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=lncRNA_filt_js
#SBATCH --output=logs/lncRNA_filt_js.%j.out
#SBATCH --time=001:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# ── Job-specific settings ────────────────────────────────────────────────────
SCREEN_NAME=lncRNA_filt_js
GENE_REF=/scratch/alpine/wist9668/reference/annotation/lnc_RNA/lncRNA_introns_no_PC_overlap_CM.bed
REF_ID=lncRNA_noPC
# ─────────────────────────────────────────────────────────────────────────────

# Clean up previous runs
rm -rf /scratch/alpine/wist9668/tmp/${SCREEN_NAME}_output/
rm -rf /scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens/${SCREEN_NAME}_output/

# Array of input files mapped to replicate numbers
FILES=(
    "/scratch/alpine/wist9668/js_hap1_control/SRR663777.fastq.gz"
)

# Run each replicate
for REP in {1..1}; do
    echo "==============================================="
    echo "Starting replicate ${REP} at $(date)"
    echo "==============================================="

    ./analyze_sli.sh -C -R=${REP} -S=${FILES[$((REP-1))]} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}

    echo "Completed replicate ${REP} at $(date)"
    echo ""
done

echo "All replicates completed at $(date)"
