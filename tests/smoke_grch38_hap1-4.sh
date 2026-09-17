#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=atesting
#SBATCH --qos=testing
#SBATCH --job-name=smoke_grch38
#SBATCH --output=logs/smoke_grch38.%j.out
#SBATCH --time=00:50:00
#SBATCH --nodes=1
#SBATCH --ntasks=16

#####################################################################
# End-to-end smoke test of the standard GRCh38 run
#####################################################################
# Runs one replicate of batch_job_hap1-4_grch38.sh's configuration
# against the 1 MB subset FASTQ, so the whole chain - env activation,
# preflight, trim, bowtie, intersectBed, binomial_test.py, compress,
# move, archive - is exercised in a few minutes instead of a few hours.
#
# atesting/testing is capped at 16 CPUs, 1 hour and one job at a time.
#
# Submit:  sbatch tests/smoke_grch38_hap1-4.sh
# Clean up afterwards (the test writes a real screen directory):
#   rm -rf /scratch/alpine/wist9668/haploid_scratch/screens/smoke_grch38_hap1-4_output
#   rm -rf /pl/active/ShenLab_PL/haploid_genetics/outputs/smoke_grch38_hap1-4_output
#####################################################################

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

SCREEN_NAME=smoke_grch38_hap1-4
GENE_REF=/pl/active/ShenLab_PL/haploid_genetics/reference/annotation/GRCh38_introns_CM_mapped.bed
REF_ID=GRCh38-NCBI_RefSeq
SEQ=/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047161_subset.fastq.gz

TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)

echo "=== Environment ==="
echo "conda env : ${CONDA_DEFAULT_ENV}"
echo "python    : $(python3 --version 2>&1)"
echo "bowtie    : $(awk -F '\t' '$1 ~ /location_bowtie/ {print $2}' alignment_settings.conf)"
echo "intersectBed : $(command -v intersectBed || echo 'NOT ON PATH')"
echo "bedtools  : $(bedtools --version 2>&1)"
echo "tmp_dir   : ${TMP_DIR}"
echo "final_dir : ${FINAL_DIR}"

mkdir -p "$TMP_DIR" "$FINAL_DIR"
rm -rf "${TMP_DIR}${SCREEN_NAME}_output/" "${FINAL_DIR}${SCREEN_NAME}_output/"
rm -rf "/pl/active/ShenLab_PL/haploid_genetics/outputs/${SCREEN_NAME}_output"

echo "=== Running one replicate on the subset FASTQ ==="
./analyze_sli.sh -C -R=1 -S=${SEQ} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}

echo "=== Output tree ==="
find "${FINAL_DIR}${SCREEN_NAME}_output" -type f -printf '%10s  %p\n' | sort -k2

echo "=== sense_vs_antisense_counts.csv head ==="
head -5 "${FINAL_DIR}${SCREEN_NAME}_output/replicate_1/${REF_ID}/sense_vs_antisense_counts.csv"
echo "=== row count ==="
wc -l < "${FINAL_DIR}${SCREEN_NAME}_output/replicate_1/${REF_ID}/sense_vs_antisense_counts.csv"

echo "=== Testing the archive step ==="
bash sub/archive_screen.sh "${SCREEN_NAME}" "${REF_ID}" "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out"

echo "=== Smoke test finished at $(date) ==="
