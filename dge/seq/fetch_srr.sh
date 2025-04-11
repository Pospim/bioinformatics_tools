#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <project_dir>"
    exit 1
fi

# Assign project directory path to variable
PROJECT_DIR="$1"

# Define subdirectories for the project
SRA_DIR="$PROJECT_DIR/sra_dir"
RAW_DATA_DIR="$PROJECT_DIR/raw_data"
QC_DIR="$PROJECT_DIR/raw_data_qc"
FINAL_DIR="$PROJECT_DIR/finfetch"
THREADS=72

# Create the necessary directories if they don't exist
mkdir -p "$SRA_DIR" "$RAW_DATA_DIR" "$QC_DIR" "$FINAL_DIR"

# Read SRR IDs from file (assumed one SRR per line)
srr_ids=($(cat "$PROJECT_DIR/srr_list.txt"))

# Set the maximum number of parallel jobs
MAX_JOBS=8

#function resubmit {
#    if [ ! -f "$RAWDIR/finished.dat" ]; then
#       sbatch --export=ALL,PRODIR="$PRODIR",SRADIR="$SRADIR",RAWDIR="$RAWDIR",QCDIR="$QCDIR",FINDIR="$FINDIR",job_restarted=1 $PRODIR/fetch_srr_list.sh $PRODIR
#    fi
#}
#trap resubmit EXIT

# Function to download and verify an SRR file
download_and_verify() {
    local srr_id=$1
    echo "Attempting to download $srr_id ..."

    # Download SRR file
    if ! prefetch "$srr_id" -O "$SRA_DIR" --max-size 50g; then
        echo "Error: Failed to download $srr_id"
        return 1
    fi
    
    # Convert to FASTQ
    if ! fastq-dump --split-files --gzip --outdir "$RAW_DATA_DIR" "$SRA_DIR/$srr_id/$srr_id.sra"; then
        echo "Error: Failed to convert $srr_id to FASTQ"
        return 1
    fi
    
     # Run FastQC and check for paired-end data
    if [[ -e "$RAW_DATA_DIR/${srr_id}_1.fastq.gz" && -s "$RAW_DATA_DIR/${srr_id}_1.fastq.gz" ]]; then
        fastqc -o "$QC_DIR" -t "$THREADS" "$RAW_DATA_DIR/${srr_id}_1.fastq.gz"
        echo "QC passed for ${srr_id}_1.fastq.gz"
        touch "$FINAL_DIR/${srr_id}.dat"
    fi
    
    if [[ -e "$RAW_DATA_DIR/${srr_id}_2.fastq.gz" && -s "$RAW_DATA_DIR/${srr_id}_2.fastq.gz" ]]; then
        fastqc -o "$QC_DIR" -t "$THREADS" "$RAW_DATA_DIR/${srr_id}_2.fastq.gz"
        echo "QC passed for ${srr_id}_2.fastq.gz"
    fi
    
    echo "------------------------------------------------------------------------"
}

# Parallelize the download process
for srr_id in "${srr_ids[@]}"; do
    if [[ ! -e "$FINAL_DIR/${srr_id}.dat" ]]; then
    	rm -rf "$SRA_DIR/${srr}"* "$RAW_DATA_DIR/${srr}"* "$QC_DIR/${srr_id}"*
    	
         # Run download and verification in parallel
        download_and_verify "$srr_id" 
        
         # Control the number of parallel jobs
        if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
            wait -n  # Wait for at least one job to finish before starting another
        fi
    fi
done
wait  # Wait for all background jobs to complete

# Mark the download process as complete
echo "Download complete."
touch "$RAW_DATA_DIR/finished.dat"

# Submit trimming job based on paired-end detection
cd "$PROJECT_DIR"

if ls "$RAW_DATA_DIR"/*_2.fastq.gz 1> /dev/null 2>&1; then
    cp "$HOME/scripts/trim_paired.sh" .
    sbatch --time=4:00:00 --mem=128G --partition=b64_any trim_paired.sh "$PROJECT_DIR"
else
    cp "$HOME/scripts/trim_single.sh" .
    sbatch --time=4:00:00 --mem=128G --partition=b64_any trim_single.sh "$PROJECT_DIR"
fi

