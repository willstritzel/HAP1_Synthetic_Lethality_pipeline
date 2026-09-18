#!/bin/bash
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=catalogue_hap1-4
#SBATCH --output=logs/catalogue_hap1-4.%j.out
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

#####################################################################
# HAP1 gene-trap screen: the four unselected HAP1 replicates
# (Blomen et al. 2015, SRR2047158-161) against the introns of the
# combined ORF screen catalogue.
#
# Feature set: combined_screen_catalogue_introns.bed - 1,966 introns
# from the 1,618 multi-exon ORFs in combined_screen_catalogue_Hsu.xlsx
# (4,845 ORFs total; single-block ORFs have no intron and contribute no
# rows). An intron is the gap between consecutive exon blocks of one
# ORF, taken from the 'starts (0-based)' / 'ends (0-based)' columns and
# written as chromStart = previous block end, chromEnd = next block
# start. Chromosome names are GenBank CM accessions from the GRCh38
# assembly report, matching the bowtie index built from
# GCA_000001405.15_GRCh38_genomic.fna.
#
# Same structure as batch_job_hap1-4_grch38.sh; only the feature set
# differs, so the sense/antisense readout here is directly comparable
# to the GRCh38_hap1-4 reference run.
#
# -C (control-creation mode) is correct here: hap1-4 IS the unselected
# control dataset, so there is nothing to normalize against and
# normalize.py / test_replicate_vs_control.py are rightly skipped.
# See the F07 note in the project docs.
#####################################################################

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# -- Job-specific settings ---------------------------------------------------
SCREEN_NAME=lncMP_combined_screen_catalogue_hap1-4
GENE_REF=/pl/active/ShenLab_PL/haploid_genetics/reference/annotation/combined_screen_catalogue_introns.bed
REF_ID=combined_screen_catalogue
# ----------------------------------------------------------------------------

# Output locations come from the settings file, so a change there does not
# have to be mirrored into every wrapper by hand.
TMP_DIR=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)
FINAL_DIR=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)

# Array of input files mapped to replicate numbers
FILES=(
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047158.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047159.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047160.fastq.gz"
    "/pl/active/ShenLab_PL/haploid_genetics/raw_reads/hap1-4/SRR2047161.fastq.gz"
)

# -- Preflight ---------------------------------------------------------------
# Fail before the first alignment rather than four hours in. The 2026 scratch
# purge left every wrapper pointing at annotation BEDs and FASTQs that no
# longer existed (F04, F05).
echo "Preflight checks:"
MISSING=0
if [[ -s $GENE_REF ]]; then
    echo "  gene_ref OK   ($(wc -l < "$GENE_REF") features) $GENE_REF"
else
    echo "  gene_ref MISSING OR EMPTY: $GENE_REF"; MISSING=1
fi
# analyze_sli.sh reads $9 as feature name and $11 as strand, so the BED must
# have exactly six columns on every line.
if [[ -s $GENE_REF ]]; then
    BADCOLS=$(awk -F '\t' 'NF!=6' "$GENE_REF" | wc -l)
    BADCOORD=$(awk -F '\t' '$2<0 || $3<=$2' "$GENE_REF" | wc -l)
    BADSTRAND=$(awk -F '\t' '$6!="+" && $6!="-"' "$GENE_REF" | wc -l)
    echo "  gene_ref format: $BADCOLS non-6-column, $BADCOORD bad-coordinate, $BADSTRAND bad-strand lines"
    [[ $BADCOLS -eq 0 && $BADCOORD -eq 0 && $BADSTRAND -eq 0 ]] || MISSING=1
    # Every BED chromosome must exist in the bowtie index, or it silently
    # intersects nothing - the CM/NC naming trap.
    FAI=/pl/active/ShenLab_PL/haploid_genetics/reference/GCA_000001405.15_GRCh38_genomic.fna.fai
    ORPHAN=$(comm -23 <(cut -f1 "$GENE_REF" | sort -u) <(cut -f1 "$FAI" | sort -u))
    if [[ -n $ORPHAN ]]; then
        echo "  gene_ref chromosomes NOT in bowtie index:"; echo "$ORPHAN"; MISSING=1
    else
        echo "  gene_ref chromosomes OK (all present in $(basename "$FAI"))"
    fi
fi
for f in "${FILES[@]}"; do
    if [[ -s $f ]]; then echo "  fastq OK      $f"
    else echo "  fastq MISSING: $f"; MISSING=1; fi
done
if [[ $MISSING -ne 0 ]]; then
    echo "Preflight failed - nothing submitted to bowtie. Exiting."; exit 1
fi
echo "  tmp_dir       $TMP_DIR"
echo "  final_dir     $FINAL_DIR"
mkdir -p "$TMP_DIR" "$FINAL_DIR"

# Clean up previous runs
rm -rf "${TMP_DIR}${SCREEN_NAME}_output/"
rm -rf "${FINAL_DIR}${SCREEN_NAME}_output/"

# Run each replicate
for REP in {1..4}; do
    echo "==============================================="
    echo "Starting replicate ${REP} at $(date)"
    echo "==============================================="

    ./analyze_sli.sh -C -R=${REP} -S=${FILES[$((REP-1))]} --name=${SCREEN_NAME} --gene-ref=${GENE_REF} --ref-id=${REF_ID}

    echo "Completed replicate ${REP} at $(date)"
    echo ""
done

echo "All replicates completed at $(date)"

# -- Archive to PetaLibrary --------------------------------------------------
# final_dir is on scratch and will be purged. Do not remove this step.
bash sub/archive_screen.sh "${SCREEN_NAME}" "${REF_ID}" "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out"
