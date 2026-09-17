#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wist9668@colorado.edu
#
# --array, --job-name and --output are set by the launcher from the CONFIG
# block below, so they are deliberately not directives here.

#####################################################################
# Screen template: one replicate per array task.
#####################################################################
# To run a new screen, copy this file, edit the CONFIG block, submit:
#
#   cp batch_job_template.sh batch_job_<annotation>_<reads>.sh
#   bash sub/submit_screen.sh batch_job_<annotation>_<reads>.sh
#
# Nothing outside CONFIG needs changing to point at a different BED or a
# different read set. The launcher reads CONFIG, sizes the array to the
# number of FASTQs, validates the annotation before anything queues,
# clears the previous run once, and chains the PetaLibrary archive
# behind the array with afterok.
#
# Do NOT sbatch this file directly. It has no cleanup step, because a
# per-task `rm -rf` of the screen directory would delete the other
# tasks' work, and no archive step, because that must wait for every
# replicate.
#
# Walltime is per replicate: measured 15-51 min over nine replicates, so
# 2 h is ~2.4x the worst case. A screen finishes in one replicate's time
# rather than the sum, at the cost of 16 CPUs per replicate and that many
# concurrent reads off the shared lab PetaLibrary.
#####################################################################

# >>> CONFIG ---------------------------------------------------------
# SCREEN_NAME  <annotation>_<reads>; output lands in <SCREEN_NAME>_output
SCREEN_NAME=GRCh38_hap1-4

# GENE_REF     BED6 feature file. Chromosome names must match the bowtie
#              index (GenBank CM accessions) or it intersects nothing.
GENE_REF=/pl/active/ShenLab_PL/haploid_genetics/reference/annotation/GRCh38_introns_CM_mapped.bed

# REF_ID       output subdirectory created inside each replicate
REF_ID=GRCh38-NCBI_RefSeq

# FILES        one FASTQ per replicate, in replicate order. The array is
#              sized to this list, so add or remove entries freely.
FILES=(
    /pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047158.fastq.gz
    /pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047159.fastq.gz
    /pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047160.fastq.gz
    /pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047161.fastq.gz
)
# <<< CONFIG ---------------------------------------------------------

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1
cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

REP=${SLURM_ARRAY_TASK_ID:?submit this through sub/submit_screen.sh, not sbatch}
SEQ=${FILES[$((REP-1))]}

TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)
mkdir -p "$TMP_DIR" "$FINAL_DIR"        # mkdir -p is concurrency-safe

echo "${SCREEN_NAME} replicate ${REP} (task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID}) on $(hostname) at $(date)"

# The launcher validated the annotation in full; re-check existence only,
# in case a filesystem was purged between submission and start.
for f in "$GENE_REF" "$SEQ"; do
    [[ -s $f ]] || { echo "MISSING: $f"; exit 1; }
done

# -C is control-creation mode. It is correct for an unselected control
# dataset, which is what every screen here currently is; see the F07 note
# in the project docs before changing it.
./analyze_sli.sh -C -R=${REP} -S="${SEQ}" --name="${SCREEN_NAME}" \
    --gene-ref="${GENE_REF}" --ref-id="${REF_ID}"
RC=$?

# Exit non-zero unless the replicate actually landed, so the afterok
# archive job cannot archive a short screen. Without this a replicate that
# dies mid-run leaves analyze_sli.sh's status unchecked and the loss
# unnoticed until someone counts the archived directories.
if [[ $RC -ne 0 || ! -d ${FINAL_DIR}${SCREEN_NAME}_output/replicate_${REP} ]]; then
    echo "FAILED: replicate ${REP} did not reach final_dir (analyze_sli.sh exit ${RC}) at $(date)"
    exit 1
fi
echo "Completed replicate ${REP} at $(date)"
