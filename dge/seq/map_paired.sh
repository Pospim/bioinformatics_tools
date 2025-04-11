#!/bin/bash
#SBATCH --mem=800G
#SBATCH --partition=b64_1024_extra
#SBATCH --time=0-168:00:00
#SBATCH -N 1

# Check if exactly one argument (project directory) is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <project_dir>"
    exit 1
fi

# Define directories
PROJECT_DIR=$1
TRIM_DATA_DIR="$PROJECT_DIR/raw_data_trim"
MAP_DATA_DIR="$PROJECT_DIR/star_mapped"
GENOME_DIR="$HOME/reference"
FINAL_DIR="$PROJECT_DIR/finmap"
THREADS=72

# Create directories if they don't exist
mkdir -p "$MAP_DATA_DIR" "$FINAL_DIR"

echo "--Mapping references--"

# Function to perform the STAR mapping
map_sample() {
    local R1="$1"
    local base_name=$(basename "$R1" _1.paired.fastq.gz)
    local output_prefix="${MAP_DATA_DIR}/${base_name}_"
    local output_sam="${output_prefix}Aligned.out.sam"
    local R2="${TRIM_DATA_DIR}/${base_name}_2.paired.fastq.gz"

    # Check if the mapping is already done
    if [[ -e "$FINAL_DIR/${base_name}.dat" ]]; then
        echo "Skipping $base_name, already mapped."
        return 0
    fi

    echo "-- Mapping sample $base_name --"
    date

    # Run STAR for paired-end mapping
    STAR \
    --runThreadN $THREADS \
    --genomeDir "$GENOME_DIR/star_index" \
    --outFileNamePrefix "$output_prefix" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --outSAMunmapped Within \
    --outSAMattributes Standard

    # Check if STAR completed successfully
    if [ $? -ne 0 ]; then
        echo "Mapping failed for sample $base_name"
        return 1
    else
        echo "Mapping completed for sample $base_name"
	echo "-------------------------------------------------------"
        # Mark the sample as successfully mapped
        touch "$FINAL_DIR/${base_name}.dat"
    fi
}

# Loop through each paired-end read file and start mapping in parallel
for R1 in "$TRIM_DATA_DIR"/*_1.paired.fastq.gz; do

    map_sample "$R1" 
done

# Mapping completed
echo "All samples have been successfully mapped."

# Submit the counting script after mapping is done
cp "$HOME/scripts/count_paired.sh" "$PROJECT_DIR"
sbatch "$PROJECT_DIR/count_paired.sh" "$PROJECT_DIR"

