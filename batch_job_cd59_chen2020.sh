#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=js-hap1-CD59-chen2020
#SBATCH --output=logs/js-hap1-CD59-chen2020.%j.out
#SBATCH --time=002:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

SCREEN_NAME=JS_HAP1_CD59_chen2020
GENE_REF=/scratch/alpine/wist9668/reference/annotation/chen2020/chen2020_intron_CM-mapped.bed
REF_ID=chen2020
SEQFILE=/pl/active/ShenLab_PL/haploid_genetics/raw_reads/CD59/SRR2016930.fastq.gz

echo "Starting CD59 chen2020 alignment at $(date)"
./analyze_sli.sh -C -R=1 -S=${SEQFILE} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}
echo "Completed at $(date)"
