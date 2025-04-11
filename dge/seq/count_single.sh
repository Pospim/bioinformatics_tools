#!/bin/bash
#SBATCH --mem=800G
#SBATCH --partition=b64_1024_long
#SBATCH --time=12:00:00
#SBATCH -N 1  

# Check if exactly one argument (project directory) is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <project_dir>"
    exit 1
fi

# Define the project directory and necessary paths
PROJECT_DIR="$1"
GENOME_DIR="$HOME/reference"
MAPPING_DIR="$PROJECT_DIR/star_mapped" 
FINAL_DIR="$PROJECT_DIR/fincount"
GTF_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.113.gtf"
CORES=64

# Create necessary directories if they don't exist
mkdir -p "$FINAL_DIR"

# Transform SAM files to sorted BAM files
for SAM in "$MAPPING_DIR"/*.sam; do
    date 

    # Define BAM file name by replacing the .sam extension with .bam
    BAM="${SAM%.sam}.bam"
    srr_id=$(basename "$SAM" | sed 's/_Aligned.*//') # Extract SRR ID from the SAM file name

    # Check if the output marker file already exists to skip processing
    if [[ -e "$FINAL_DIR/${srr_id}.dat" ]]; then
        echo "Skipping $srr_id as it already exists"
        echo "---------------------------------------------------"
        continue
    else
        echo "Transforming $srr_id to sorted BAM"
        # Convert SAM to BAM and sort it using samtools
        samtools view -Sb "$SAM" | samtools sort -@ "$CORES" -m 10g > "$BAM"
    fi

    # Check if the transformation was successful
    if [ $? -ne 0 ]; then
        echo "Error: Transformation failed for $SAM"
        exit 1
    fi

    # Create a marker file indicating successful transformation
    touch "$FINAL_DIR/${srr_id}.dat"
    echo "-----------------------------------------------------------"
done

date
echo "-- All SAM files transformed to sorted BAM --"

# Prepare a read matrix table using featureCounts for single-end data
featureCounts -T "$CORES" -M -a "$GTF_FILE" -o "$PROJECT_DIR/featurecounts_all.txt" "$MAPPING_DIR"/*.bam 

if [ $? -ne 0 ]; then
    echo "Error: featureCounts failed"
    exit 1
fi

date
echo "-- Read matrix table prepared using featureCounts --"

