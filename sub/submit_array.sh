#!/bin/bash
#####################################################################
# Launch a parallel screen: clean once, submit the array, chain the
# archive job behind it.
#####################################################################
# This is a plain shell script, NOT an sbatch script - run it on a
# login node (the DTN's Slurm client cannot reach the controller).
#
# It exists because the two things that must happen exactly once per
# screen cannot live inside an array task:
#   - removing the previous run's output (a per-task rm -rf would
#     delete the other tasks' work)
#   - archiving to the PetaLibrary (must wait for every replicate)
#
# Usage:
#   bash sub/submit_array.sh <array_wrapper.sh> <SCREEN_NAME> <REF_ID>
#
# Example:
#   bash sub/submit_array.sh batch_job_hap1-4_grch38_array.sh \
#        ControlData-HAP1 GRCh38-NCBI_RefSeq
#####################################################################

set -euo pipefail

WRAPPER=${1:?usage: submit_array.sh <array_wrapper.sh> <SCREEN_NAME> <REF_ID>}
SCREEN_NAME=${2:?usage: submit_array.sh <array_wrapper.sh> <SCREEN_NAME> <REF_ID>}
REF_ID=${3:?usage: submit_array.sh <array_wrapper.sh> <SCREEN_NAME> <REF_ID>}

[[ -f $WRAPPER ]] || { echo "No such wrapper: $WRAPPER"; exit 1; }
[[ -f alignment_settings.conf ]] || { echo "Run this from the pipeline directory."; exit 1; }

if ! sinfo -h -o '%P' >/dev/null 2>&1; then
    echo "Slurm is not reachable from this host. Run this on a login node, not the DTN."
    exit 1
fi

TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)

echo "Screen     : ${SCREEN_NAME}"
echo "Reference  : ${REF_ID}"
echo "tmp_dir    : ${TMP_DIR}"
echo "final_dir  : ${FINAL_DIR}"

# Clean up the previous run - once, before anything is submitted.
rm -rf "${TMP_DIR}${SCREEN_NAME}_output/"
rm -rf "${FINAL_DIR}${SCREEN_NAME}_output/"
mkdir -p "$TMP_DIR" "$FINAL_DIR" logs

# sbatch prints a partition-retirement banner, so pull the id out by pattern
# rather than assuming the id is the only thing on stdout.
AID=$(sbatch "$WRAPPER" | grep -oE '[0-9]{5,}' | tail -1)
[[ -n $AID ]] || { echo "Could not read the array job id from sbatch."; exit 1; }
echo "Array job  : ${AID}"

# afterok fires only if EVERY array task exits 0, so a failed replicate leaves
# the archive job pending rather than archiving a partial screen.
BID=$(sbatch --dependency=afterok:"${AID}" \
             --job-name="archive_${SCREEN_NAME}" \
             --output="logs/archive_${SCREEN_NAME}.%j.out" \
             sub/archive_job.sh "${SCREEN_NAME}" "${REF_ID}" "${AID}" \
      | grep -oE '[0-9]{5,}' | tail -1)
echo "Archive job: ${BID} (waits for ${AID} to finish clean)"
echo
echo "Watch it with:  squeue -u \$USER"
echo "If a replicate fails, the archive job stays PENDING with reason"
echo "DependencyNeverSatisfied; cancel it with:  scancel ${BID}"
