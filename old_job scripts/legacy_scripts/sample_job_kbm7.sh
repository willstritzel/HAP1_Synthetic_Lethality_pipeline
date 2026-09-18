#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=kbm7-controls
#SBATCH --output=kbm7-controls.%j.out
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline


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
    
    ./analyze_sli.sh -C -R=${REP} -S=${FILES[$((REP-1))]} --name=ControlData-KBM7
    
    echo "Completed replicate ${REP} at $(date)"
    echo ""
done

echo "All replicates completed at $(date)"
