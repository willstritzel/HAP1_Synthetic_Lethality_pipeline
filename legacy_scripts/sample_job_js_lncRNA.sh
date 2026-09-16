#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=js-lncRNA-control
#SBATCH --output=js-lncRNA-control.%j.out
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline
# Clean up previous runs
rm -rf /scratch/alpine/wist9668/tmp/JS_ControlData-HAP1-lncRNA_output/
rm -rf /scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens/JS_ControlData-HAP1-lncRNA_output/

echo "==============================================="
echo "Starting JS lncRNA control analysis at $(date)"
echo "==============================================="

./analyze_sli.sh -C -R=1 -S=/scratch/alpine/wist9668/js_hap1_control/SRR663777.fastq.gz --name=JS_ControlData-HAP1-lncRNA

echo "==============================================="
echo "Completed at $(date)"
echo "==============================================="
