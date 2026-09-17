#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=smoke_arr
#SBATCH --output=logs/smoke_arr.%A_%a.out
#SBATCH --array=1-2
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=16

#####################################################################
# Concurrency test for the analyze_sli.sh parallel-safety fixes.
#####################################################################
# Runs TWO replicates of one screen at the same time on the 1 MB subset
# FASTQ. Under the pre-2026-09-16 code this could not work: whichever
# task started second hit the screen-level tmp guard and exited, and
# whichever finished first deleted the shared tmp parent out from under
# the other. Both tasks completing with their own output directory is
# the pass condition.
#
# Deliberately on acpu rather than atesting: the testing QOS allows one
# job at a time, which is the opposite of what needs testing here.
#
# Submit:  sbatch tests/smoke_array_grch38.sh
#####################################################################

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

SCREEN_NAME=smoke_arr_grch38
GENE_REF=/pl/active/ShenLab_PL/haploid_genetics/reference/annotation/GRCh38_introns_CM_mapped.bed
REF_ID=GRCh38-NCBI_RefSeq
SEQ=/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047161_subset.fastq.gz

REP=${SLURM_ARRAY_TASK_ID}

TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)
mkdir -p "$TMP_DIR" "$FINAL_DIR"

echo "task ${SLURM_ARRAY_TASK_ID} / job ${SLURM_ARRAY_JOB_ID}: replicate ${REP} on $(hostname) at $(date)"
./analyze_sli.sh -C -R=${REP} -S=${SEQ} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}
echo "task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
