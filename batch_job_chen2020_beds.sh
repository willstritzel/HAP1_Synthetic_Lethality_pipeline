#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=chen2020-beds
#SBATCH --output=logs/chen2020-beds.%j.out
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

echo "Starting Chen et al. 2020 BED creation at $(date)"
conda run -n hap1 python3 sub/make_chen2020_beds_CM_mapped.py
echo "Completed at $(date)"
