#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --job-name=tiled_2kb_hap1-4
#SBATCH --output=logs/tiled_2kb_hap1-4.%j.out
#SBATCH --time=002:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wist9668@colorado.edu

module purge
module load miniforge
mamba activate hap1

cd /projects/wist9668/HAP1_Synthetic_Lethality_pipeline

# Set gene_ref and ref_id for 2kb tiled windows (1kb lead offset)
sed -i 's|^gene_ref.*|gene_ref\t/scratch/alpine/wist9668/reference/annotation/tiled_windows/genome_2kb_windows_1kb_offset.bed|' alignment_settings.conf
sed -i 's|^ref_id.*|ref_id\ttiled_2kb_hap1-4|' alignment_settings.conf

# Clean up previous runs
rm -rf /scratch/alpine/wist9668/tmp/tiled_2kb_hap1-4_output/
rm -rf /scratch/alpine/wist9668/data/analyzed_data/synthetic_lethal_screens/tiled_2kb_hap1-4_output


# Array of input files mapped to replicate numbers
FILES=(
    "/scratch/alpine/wist9668/hap1-4/SRR2047158.fastq.gz"
    "/scratch/alpine/wist9668/hap1-4/SRR2047159.fastq.gz"
    "/scratch/alpine/wist9668/hap1-4/SRR2047160.fastq.gz"
    "/scratch/alpine/wist9668/hap1-4/SRR2047161.fastq.gz"
)

# Run each replicate
for REP in {1..4}; do
    echo "==============================================="
    echo "Starting replicate ${REP} at $(date)"
    echo "==============================================="

    ./analyze_sli.sh -C -R=${REP} -S=${FILES[$((REP-1))]} --name=tiled_2kb_hap1-4

    echo "Completed replicate ${REP} at $(date)"
    echo ""
done

echo "All replicates completed at $(date)"
