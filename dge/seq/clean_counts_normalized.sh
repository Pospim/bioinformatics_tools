#!/bin/bash

# Check if the number of arguments is 1
if [ $# -ne 1 ]; then
    echo "Usage: $0 <project_dir>"
    exit 1
fi

PRODIR=$1
PROJECT=$(basename "$PRODIR")

if ! cd "$PRODIR"; then
    echo "Error: invalid directory"
    exit 1
fi

sed '1d' featurecounts_all.txt > featurecounts.txt 

if [[ "$PROJECT" =~ ^(GSE[0-9]+) ]]; then
    GSE="${BASH_REMATCH[1]}"
    sed -e "1s|/home2/pospisil/projects/$GSE/star_mapped/||g" -e '1s/_Aligned\.out\.bam//g' featurecounts.txt > cleaned_counts.txt
else
    sed -E 's|([^[:space:]]*/)?(SRR[0-9]+)_Aligned[^[:space:]]*|\2|g' featurecounts.txt > cleaned_counts.txt
fi  # <-- Missing fi to close the if-else block

cut --complement -f 2-5 cleaned_counts.txt > counts.txt

python3 /home/pospim/Desktop/work/bioinformatics/tools/dge/normalize_counts.py -i counts.txt

cut -f1,3- counts_normalized.txt | sponge counts_normalized.csv

echo "Counts extracted successfully"

