#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=chen2020-beds
#SBATCH --output=logs/chen2020-beds.%j.out
#SBATCH --time=002:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

echo "Starting Chen et al. 2020 BED creation at $(date)"
conda run -n hap1 python3 sub/make_chen2020_beds_CM_mapped.py
echo "Completed at $(date)"
