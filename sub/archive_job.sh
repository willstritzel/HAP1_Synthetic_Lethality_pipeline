#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

#####################################################################
# Slurm wrapper around sub/archive_screen.sh, for the parallel path.
#####################################################################
# The serial wrappers call archive_screen.sh inline at the end of the
# job. An array cannot do that - no single task knows the others have
# finished - so the launcher submits this with --dependency=afterok on
# the array.
#
# One task is right: this is pure I/O. The hour of walltime is for the
# PetaLibrary mount, which has been observed to slow to a crawl; the
# copy itself is ~6 GB.
#
# Usage (normally via sub/submit_array.sh):
#   sbatch sub/archive_job.sh <SCREEN_NAME> <REF_ID> [ARRAY_JOB_ID]
#####################################################################

SCREEN_NAME=${1:?usage: sbatch sub/archive_job.sh <SCREEN_NAME> <REF_ID> [ARRAY_JOB_ID]}
REF_ID=${2:?usage: sbatch sub/archive_job.sh <SCREEN_NAME> <REF_ID> [ARRAY_JOB_ID]}
ARRAY_JOB_ID=${3:-}

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

bash sub/archive_screen.sh "${SCREEN_NAME}" "${REF_ID}"

# Archive the array's own Slurm logs alongside the data. The output pattern is
# %A_%a, so there is one file per replicate; match on the array job id rather
# than a bare wildcard, which would drag in other screens' logs.
DEST=/pl/active/ShenLab_PL/haploid_genetics/outputs/${SCREEN_NAME}_output
if [[ -n $ARRAY_JOB_ID ]]; then
    for f in logs/*."${ARRAY_JOB_ID}"_*.out; do
        [[ -f $f ]] && cp -p "$f" "${DEST}/"
    done
fi

echo "=== Archived replicates ==="
find "${DEST}" -name 'sense_vs_antisense_counts.csv' | sort
