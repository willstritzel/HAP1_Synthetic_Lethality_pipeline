#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=JOBNAME
#SBATCH --output=logs/JOBNAME.%j.out
#SBATCH --time=002:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# ── Job-specific settings ────────────────────────────────────────────────────
SCREEN_NAME=SCREENNAME
GENE_REF=/path/to/gene_ref.bed
REF_ID=ref_id_label
# ─────────────────────────────────────────────────────────────────────────────

# Clean up previous runs
rm -rf /scratch/alpine/wist9668/tmp/${SCREEN_NAME}_output/
rm -rf /scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens/${SCREEN_NAME}_output/

# Array of input files mapped to replicate numbers
FILES=(
    "/scratch/alpine/wist9668/DATASET/SRR_REP1.fastq.gz"
    "/scratch/alpine/wist9668/DATASET/SRR_REP2.fastq.gz"
    "/scratch/alpine/wist9668/DATASET/SRR_REP3.fastq.gz"
    "/scratch/alpine/wist9668/DATASET/SRR_REP4.fastq.gz"
)

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
