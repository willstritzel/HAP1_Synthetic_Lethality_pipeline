#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16

#####################################################################
# Re-run ONE replicate of an existing screen.
#####################################################################
# The batch_job_* wrappers begin by deleting the whole screen
# directory, so they cannot be used to fill in a single missing
# replicate - they would destroy the ones that succeeded.
#
# A replicate goes missing when its analyze_sli.sh dies mid-run: the
# wrapper's loop prints "Completed replicate N" without checking the
# exit status, and archive_screen.sh then archives whatever it finds.
# Compare the replicate_* directories in final_dir against the number
# of FASTQs to spot it.
#
# This script deletes only the named replicate's leftover working
# directory, then runs that one replicate. It does NOT archive; run
# sub/archive_screen.sh afterwards, which re-copies the whole screen.
#
# Usage:
#   sbatch sub/rerun_replicate.sh <SCREEN_NAME> <REP> <FASTQ> <GENE_REF> <REF_ID>
#####################################################################

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

SCREEN_NAME=${1:?usage: sbatch sub/rerun_replicate.sh <SCREEN_NAME> <REP> <FASTQ> <GENE_REF> <REF_ID>}
REP=${2:?missing REP}
SEQ=${3:?missing FASTQ}
GENE_REF=${4:?missing GENE_REF}
REF_ID=${5:?missing REF_ID}

TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)

echo "Re-running ${SCREEN_NAME} replicate ${REP} on $(hostname) at $(date)"
echo "  fastq    : ${SEQ}"
echo "  gene_ref : ${GENE_REF}"
echo "  ref_id   : ${REF_ID}"
echo "  final_dir: ${FINAL_DIR}"

MISSING=0
[[ -s $SEQ ]]      || { echo "  fastq MISSING: $SEQ"; MISSING=1; }
[[ -s $GENE_REF ]] || { echo "  gene_ref MISSING: $GENE_REF"; MISSING=1; }
if [[ -d ${FINAL_DIR}${SCREEN_NAME}_output/replicate_${REP} ]]; then
    echo "  replicate_${REP} ALREADY PRESENT in final_dir - refusing to overwrite."
    MISSING=1
fi
[[ $MISSING -eq 0 ]] || { echo "Preflight failed. Exiting."; exit 1; }

# Only this replicate's leftovers, never the screen directory.
rm -rf "${TMP_DIR}${SCREEN_NAME}_output/replicate_${REP}"
mkdir -p "$TMP_DIR" "$FINAL_DIR"

./analyze_sli.sh -C -R=${REP} -S=${SEQ} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}
RC=$?

echo "analyze_sli.sh exit status: ${RC} at $(date)"
if [[ ! -d ${FINAL_DIR}${SCREEN_NAME}_output/replicate_${REP} ]]; then
    echo "FAILED: replicate_${REP} did not reach final_dir."
    exit 1
fi
echo "Replicate ${REP} complete. Now run:"
echo "  bash sub/archive_screen.sh ${SCREEN_NAME} ${REF_ID}"
