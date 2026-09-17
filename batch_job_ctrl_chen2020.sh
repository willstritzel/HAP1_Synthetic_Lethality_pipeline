#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=js-ctrl-chen2020
#SBATCH --output=logs/js-ctrl-chen2020.%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

SCREEN_NAME=JS_ControlData_chen2020
GENE_REF=/scratch/alpine/wist9668/haploid_scratch/reference/annotation/chen2020/chen2020_intron_CM-mapped.bed
REF_ID=chen2020
SEQFILE=/scratch/alpine/wist9668/js_hap1_control/SRR663777.fastq.gz

echo "Starting js_hap1_control chen2020 alignment at $(date)"
./analyze_sli.sh -C -R=1 -S=${SEQFILE} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}
echo "Completed at $(date)"
