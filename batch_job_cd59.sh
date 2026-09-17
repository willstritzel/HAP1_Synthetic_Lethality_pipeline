#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=js-hap1-CD59
#SBATCH --output=logs/js-hap1-CD59.%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

SCREEN_NAME=JS_HAP1_CD59
GENE_REF=/scratch/alpine/wist9668/haploid_scratch/reference/annotation/GRCh38_introns_CM_mapped.bed
REF_ID=CD59
SEQFILE=/pl/active/ShenLab_PL/haploid_genetics/raw_reads/CD59/SRR2016930.fastq.gz

echo "Starting CD59 alignment at $(date)"
echo "FASTQ: ${SEQFILE}"
echo "Gene ref: ${GENE_REF}"

./analyze_sli.sh -C -R=1 -S=${SEQFILE} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}

echo "Completed at $(date)"
