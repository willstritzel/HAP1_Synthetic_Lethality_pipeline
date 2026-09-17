#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=grch38_hap1-4_arr
#SBATCH --output=logs/grch38_hap1-4_arr.%A_%a.out
#SBATCH --array=1-4
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wist9668@colorado.edu

#####################################################################
# Standard GRCh38 run, one replicate per array task.
#
# Same screen as batch_job_hap1-4_grch38.sh, which remains the simpler
# thing to reach for. This one finishes in one replicate's time
# (~40 min) instead of ~2 h 30, at the cost of 64 CPUs and four
# concurrent reads off the shared lab PetaLibrary.
#
# Walltime: the slowest replicate measured was 42m35s and the observed
# per-replicate spread on identical inputs is 15-43 min, so 1h30 is
# ~2.1x the worst case.
#
# DO NOT submit this directly - it deliberately has no cleanup step,
# because a per-task `rm -rf` of the screen directory would delete the
# other tasks' work. Submit through the launcher, which cleans once and
# chains the archive job:
#
#   bash sub/submit_array.sh batch_job_hap1-4_grch38_array.sh \
#        GRCh38_hap1-4 GRCh38-NCBI_RefSeq
#####################################################################

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# -- Job-specific settings ---------------------------------------------------
SCREEN_NAME=GRCh38_hap1-4
GENE_REF=/pl/active/ShenLab_PL/haploid_genetics/reference/annotation/GRCh38_introns_CM_mapped.bed
REF_ID=GRCh38-NCBI_RefSeq
# ----------------------------------------------------------------------------

REP=${SLURM_ARRAY_TASK_ID}

TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)

FILES=(
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047158.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047159.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047160.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047161.fastq.gz"
)
SEQ=${FILES[$((REP-1))]}

# -- Preflight ---------------------------------------------------------------
echo "Array task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID}: replicate ${REP} on $(hostname)"
MISSING=0
if [[ -s $GENE_REF ]]; then
    echo "  gene_ref OK   ($(wc -l < "$GENE_REF") features) $GENE_REF"
else
    echo "  gene_ref MISSING OR EMPTY: $GENE_REF"; MISSING=1
fi
if [[ -s $SEQ ]]; then
    echo "  fastq OK      $SEQ"
else
    echo "  fastq MISSING: $SEQ"; MISSING=1
fi
if [[ $MISSING -ne 0 ]]; then
    echo "Preflight failed - nothing submitted to bowtie. Exiting."; exit 1
fi

# mkdir -p is concurrency-safe; analyze_sli.sh needs both to exist.
mkdir -p "$TMP_DIR" "$FINAL_DIR"

echo "Starting replicate ${REP} at $(date)"
./analyze_sli.sh -C -R=${REP} -S=${SEQ} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}
echo "Completed replicate ${REP} at $(date)"
