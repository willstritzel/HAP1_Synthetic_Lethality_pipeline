#!/bin/bash
#####################################################################
# Submit a screen: validate, clear once, size the array, chain archive.
#####################################################################
# Run from the pipeline directory on a LOGIN NODE - the DTN's Slurm
# client cannot reach the controller.
#
#   bash sub/submit_screen.sh batch_job_<annotation>_<reads>.sh
#
# The wrapper's CONFIG block is the single source of truth: SCREEN_NAME,
# GENE_REF, REF_ID and the replicate count are read out of it, so they
# are never passed twice and cannot drift out of step.
#
# Three things must happen exactly once per screen and so cannot live
# inside an array task:
#   - validating the annotation (do it before anything queues)
#   - clearing the previous run (a per-task rm -rf would delete the
#     other tasks' work)
#   - archiving to the PetaLibrary (must wait for every replicate)
#####################################################################

set -euo pipefail

WRAPPER=${1:?usage: bash sub/submit_screen.sh <wrapper.sh>}
[[ -f $WRAPPER ]] || { echo "No such wrapper: $WRAPPER"; exit 1; }
[[ -f alignment_settings.conf ]] || { echo "Run this from the pipeline directory."; exit 1; }
sinfo -h -o '%P' >/dev/null 2>&1 || {
    echo "Slurm is unreachable from $(hostname) - use a login node, not the DTN."; exit 1; }

FILES=()
eval "$(sed -n '/^# >>> CONFIG/,/^# <<< CONFIG/p' "$WRAPPER")"
: "${SCREEN_NAME:?$WRAPPER has no SCREEN_NAME}"
: "${GENE_REF:?$WRAPPER has no GENE_REF}"
: "${REF_ID:?$WRAPPER has no REF_ID}"
N=${#FILES[@]}
(( N > 0 )) || { echo "$WRAPPER has an empty FILES list."; exit 1; }

TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)
REF_GENOME=$(awk -F '\t' '$1 ~ /ref_genome/ {print $2}' alignment_settings.conf)
FAI=${REF_GENOME%_bowtie}.fna.fai

printf 'Screen     : %s\nAnnotation : %s (ref_id %s)\nReplicates : %d\nfinal_dir  : %s\n\n' \
       "$SCREEN_NAME" "$(basename "$GENE_REF")" "$REF_ID" "$N" "$FINAL_DIR"

# -- Preflight -------------------------------------------------------
FAIL=0
if [[ -s $GENE_REF ]]; then
    # analyze_sli.sh reads $9 as the feature name and $11 as its strand,
    # which requires exactly six columns on every line.
    BAD=$(awk -F '\t' 'NF!=6 || $2<0 || $3<=$2 || ($6!="+" && $6!="-")' "$GENE_REF" | wc -l)
    echo "  annotation : $(wc -l < "$GENE_REF") features, $BAD malformed lines"
    (( BAD == 0 )) || FAIL=1
    # Every BED chromosome must exist in the index or it silently
    # intersects nothing - the CM/NC accession trap.
    if [[ -s $FAI ]]; then
        ORPHAN=$(comm -23 <(cut -f1 "$GENE_REF" | sort -u) <(cut -f1 "$FAI" | sort -u))
        if [[ -n $ORPHAN ]]; then
            echo "  chromosomes NOT in the bowtie index:"
            echo "$ORPHAN" | sed 's/^/               /'
            FAIL=1
        else
            echo "  chromosomes: all present in $(basename "$FAI")"
        fi
    else
        echo "  WARNING: no .fai at $FAI - chromosome naming NOT checked"
    fi
else
    echo "  annotation MISSING OR EMPTY: $GENE_REF"; FAIL=1
fi

for f in "${FILES[@]}"; do
    [[ -s $f ]] || { echo "  fastq MISSING: $f"; FAIL=1; }
done
(( FAIL == 0 )) || { echo -e "\nPreflight failed - nothing submitted."; exit 1; }
echo "  fastqs     : $N present"

# -- Clear the previous run, once, then submit -----------------------
rm -rf "${TMP_DIR}${SCREEN_NAME}_output/" "${FINAL_DIR}${SCREEN_NAME}_output/"
mkdir -p "$TMP_DIR" "$FINAL_DIR" logs

# sbatch prints a deprecation banner, so pull the id out by pattern rather
# than assuming it is the only thing on stdout.
AID=$(sbatch --array=1-"$N" \
             --job-name="$SCREEN_NAME" \
             --output="logs/${SCREEN_NAME}.%A_%a.out" \
             "$WRAPPER" | grep -oE '[0-9]{5,}' | tail -1)
[[ -n $AID ]] || { echo "Could not read the array job id from sbatch."; exit 1; }

# afterok fires only if EVERY task exits 0, and each task exits non-zero
# unless its replicate reached final_dir - so a lost replicate leaves the
# archive pending rather than archiving a short screen.
BID=$(sbatch --dependency=afterok:"$AID" \
             --job-name="archive_${SCREEN_NAME}" \
             --output="logs/archive_${SCREEN_NAME}.%j.out" \
             sub/archive_job.sh "$SCREEN_NAME" "$REF_ID" "$AID" \
      | grep -oE '[0-9]{5,}' | tail -1)

printf '\nArray job  : %s (tasks 1-%d)\nArchive job: %s (afterok:%s)\n\n' "$AID" "$N" "$BID" "$AID"
echo "Watch:   squeue -u \$USER"
echo "Logs:    logs/${SCREEN_NAME}.${AID}_<task>.out"
echo "If a replicate fails the archive stays PENDING with reason"
echo "DependencyNeverSatisfied - inspect the logs, then: scancel ${BID}"
