#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=hap1-controls
#SBATCH --output=hap1-controls.%j.out
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# Clean up previous runs
rm -rf /scratch/alpine/wist9668/tmp/ControlData-HAP1_output/
rm -rf /scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens/ControlData-HAP1_output/

# Array of input files mapped to replicate numbers
FILES=(
    "/scratch/alpine/wist9668/hap1-4/SRR2047158.fastq.gz"
    "/scratch/alpine/wist9668/hap1-4/SRR2047159.fastq.gz"
    "/scratch/alpine/wist9668/hap1-4/SRR2047160.fastq.gz"
    "/scratch/alpine/wist9668/hap1-4/SRR2047161.fastq.gz"
)

# Run each replicate
for REP in {1..4}; do
    echo "==============================================="
    echo "Starting replicate ${REP} at $(date)"
    echo "==============================================="
    
    ./analyze_sli.sh -C -R=${REP} -S=${FILES[$((REP-1))]} --name=ControlData-HAP1
    
    echo "Completed replicate ${REP} at $(date)"
    echo ""
done

echo "All replicates completed at $(date)"
