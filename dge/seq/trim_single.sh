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

# Define project directory and necessary subdirectories
PROJECT_DIR="$1"
RAW_DATA_DIR="${PROJECT_DIR}/raw_data"
TRIM_DATA_DIR="${PROJECT_DIR}/raw_data_trim"
TRIM_TOOL_DIR="$HOME/software/Trimmomatic-0.39"
QC_DIR="${PROJECT_DIR}/raw_data_trim_qc"
FINAL_DIR="${PROJECT_DIR}/fintrim"
THREADS=72

# Create the necessary directories if they don't exist
mkdir -p "$TRIM_DATA_DIR" "$QC_DIR" "$FINAL_DIR"

cd "$PROJECT_DIR"
# Function to run the trimming process
trim_reads() {
    local fastq_1="$1"
    local base_name=$(basename "$fastq_1" _1.fastq.gz)

    local output_file="$TRIM_DATA_DIR/${base_name}_1_trimmed.fastq.gz"

    # Check if output files already exist and skip if true
    if [[ -e "$FINAL_DIR/${base_name}.dat" ]]; then
        echo "Skipping trimming for $base_name as it is already processed."
        return 0
    fi

    echo "----------------------------------------------------------------"
    echo "Trimming $base_name ..."
   
    # Clean any incomplete or temporary files from previous runs
    rm -rf "$TRIM_DATA_DIR/${base_name}"* "$QC_DIR/${base_name}"*

    # Run Trimmomatic for paired-end trimming
    java -jar "$TRIM_TOOL_DIR/trimmomatic-0.39.jar" SE \
    -threads $THREADS "$fastq_1" "$output_file" \
    ILLUMINACLIP:"$TRIM_TOOL_DIR/adapters/TruSeq3-SE.fa":2:30:10:2:True \
    LEADING:3 TRAILING:3 MINLEN:36

    # Check if Trimmomatic succeeded
    if [ $? -ne 0 ]; then
        echo "Error: Trimming failed for $base_name"
        return 1
    fi
    echo "Trimming completed for $base_name."
     
    # Run FastQC for quality control
    echo "Running fastQC on $output_file"
    fastqc -o "$QC_DIR" -t $THREADS "$output_file"
    
  # Verify FastQC was successful
    if [ $? -ne 0 ]; then
        echo "Error: FastQC failed for $base_name"
        return 1
    fi

    # Mark the trimming process as complete for this sample
    touch "$FINAL_DIR/${base_name}.dat"
    echo "Trimming and QC completed for $base_name."
}

# Trim all fastq files
for fastq in "$RAW_DATA_DIR"/*_1.fastq.gz; do
   
    trim_reads "$fastq" 
    
done

# Mark the overall process as successful
echo "All reads trimmed successfully."
touch "$TRIM_DATA_DIR/finished.dat"

# Submit the mapping script based on paired-end reads
cp "$HOME/scripts/map_single.sh" "$PROJECT_DIR"
sbatch "$PROJECT_DIR/map_single.sh" "$PROJECT_DIR"

