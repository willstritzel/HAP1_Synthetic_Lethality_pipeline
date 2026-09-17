#!/bin/bash
#####################################################################
# Archive a finished screen to the lab PetaLibrary
#####################################################################
# final_dir lives on /scratch/alpine, which purges unused files after
# 90 days. That purge destroyed 22 of 32 screens in 2026 (audit F01).
# Every wrapper must therefore end by calling this script.
#
# What is copied, per replicate:
#   genome_align_log.txt          run log for the alignment stage
#   alignment_settings.conf       the settings this run actually used
#   replicate.mapped_unique.bwt   deduplicated insertions (~2.5 GB) -
#                                 the input to every sub/*.py analysis
#   <ref_id>/                     whole analysis-ready directory:
#                                 sense_vs_antisense_counts.csv,
#                                 results.csv, the gene+screen_*.txt.gz
#                                 files, annotation_analysis_log.txt and
#                                 the analytical_settings.conf copy
#
# NOT copied: replicate.fa, replicate.mapped_all.bwt,
# replicate.mapped_mm*.bwt, replicate.suppressed_mm*.bwt,
# replicate.input4enrichment.bed.gz - all regenerable from the FASTQ.
#
# Usage: bash sub/archive_screen.sh <SCREEN_NAME> <REF_ID> [SLURM_LOG]
#####################################################################

set -u

SCREEN_NAME=${1:?usage: archive_screen.sh <SCREEN_NAME> <REF_ID> [SLURM_LOG]}
REF_ID=${2:?usage: archive_screen.sh <SCREEN_NAME> <REF_ID> [SLURM_LOG]}
SLURM_LOG=${3:-}

ARCHIVE_ROOT=/pl/active/ShenLab_PL/haploid_genetics/outputs

# Read final_dir from the settings file rather than duplicating it here,
# so the archive step always tracks the configured output location.
final_dir=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)

SRC=${final_dir}${SCREEN_NAME}_output
DEST=${ARCHIVE_ROOT}/${SCREEN_NAME}_output

if [[ ! -d $SRC ]]; then
    printf "ARCHIVE FAILED: %s not found\n" "$SRC"
    exit 1
fi

printf "\n=== Archiving %s ===\n  from %s\n  to   %s\n" "$SCREEN_NAME" "$SRC" "$DEST"
mkdir -p "$DEST"

for REPDIR in "$SRC"/replicate_*; do
    [[ -d $REPDIR ]] || continue
    REP=$(basename "$REPDIR")
    mkdir -p "$DEST/$REP"

    for f in genome_align_log.txt alignment_settings.conf replicate.mapped_unique.bwt; do
        if [[ -f $REPDIR/$f ]]; then
            cp -p "$REPDIR/$f" "$DEST/$REP/$f"
        else
            printf "  WARNING: %s/%s missing\n" "$REP" "$f"
        fi
    done

    if [[ -d $REPDIR/$REF_ID ]]; then
        cp -rp "$REPDIR/$REF_ID" "$DEST/$REP/"
    else
        printf "  WARNING: %s/%s/ missing\n" "$REP" "$REF_ID"
    fi

    printf "  %s archived\n" "$REP"
done

# The Slurm log is the only record of the preflight checks and of any
# replicate that failed before writing output.
if [[ -n $SLURM_LOG && -f $SLURM_LOG ]]; then
    cp -p "$SLURM_LOG" "$DEST/$(basename "$SLURM_LOG")"
    printf "  slurm log archived: %s\n" "$(basename "$SLURM_LOG")"
fi

sync
printf "=== Archive complete: %s ===\n" "$(du -sh "$DEST" | cut -f1)"
find "$DEST" -name 'sense_vs_antisense_counts.csv' | sort
