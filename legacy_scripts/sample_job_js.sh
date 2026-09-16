#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=js-hap1-control
#SBATCH --output=js-hap1-control.%j.out
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
echo "Starting JS control analysis at $(date)"
echo "==============================================="

./analyze_sli.sh -C -R=1 -S=/scratch/alpine/wist9668/js_hap1_control/SRR663777.fastq.gz --name=JS_ControlData-HAP1

echo "==============================================="
echo "Completed at $(date)"
echo "==============================================="
