#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=js-hap1-PrP
#SBATCH --output=js-hap1-PrP.%j.out
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

echo "==============================================="
echo "Starting JS PrP analysis at $(date)"
echo "==============================================="

./analyze_sli.sh -C -R=1 -S=/scratch/alpine/wist9668/PrP/SRR2016931.fastq.gz --name=JS_HAP1_PrP

echo "==============================================="
echo "Completed at $(date)"
echo "==============================================="
