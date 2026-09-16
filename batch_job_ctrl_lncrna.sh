#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=js-ctrl-lncRNA
#SBATCH --output=logs/js-ctrl-lncRNA.%j.out
#SBATCH --time=002:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

SCREEN_NAME=JS_ControlData_lncRNA
GENE_REF=/scratch/alpine/wist9668/reference/annotation/lnc_RNA/chothani_lncRNA_intron_real_CM-mapped.bed
REF_ID=chothani_lncRNA
SEQFILE=/scratch/alpine/wist9668/js_hap1_control/SRR663777.fastq.gz

echo "Starting js_hap1_control lncRNA alignment at $(date)"
./analyze_sli.sh -C -R=1 -S=${SEQFILE} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}
echo "Completed at $(date)"
