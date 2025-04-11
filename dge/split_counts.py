import os
import csv

delim = '\t'

def parse_column_name(col_name):
    """
    Parse column name of the form:
    D[number]_[timepoint]_[condition]_[replicate number]
    For example: D9_18h_FLU_4
    Returns:
        donor (e.g. 'D9'),
        timepoint (e.g. '4h' or '18h'),
        condition (e.g. 'Mock' or 'FLU'),
        replicate (integer replicate number)
    """
    parts = col_name.split('_')
    donor = parts[0]      # e.g. D9
    timepoint = parts[1]  # e.g. 4h or 18h
    condition = parts[2]  # e.g. Mock or FLU
    replicate = int(parts[3]) if len(parts) > 3 else 1
    return donor, timepoint, condition, replicate

# Input file (replace with actual path)
input_file = '/home/pospim/Desktop/work/bioinformatics/datasets/Flu/GSE186908_/counts_named.csv'

# Derive output directory by removing the filename part
output_dir = os.path.dirname(input_file)

# Dictionary structure: data_dict[donor][timepoint] = list of (col_index)
data_dict = {}

# We need to store column info by index
col_info = {}  # col_index -> (donor, timepoint, condition, replicate)

with open(input_file, 'r', newline='') as infile:
    reader = csv.reader(infile, delimiter=delim)
    header = next(reader)
    geneid_index = 0

    # Parse column names
    for i, col_name in enumerate(header[1:], start=1):
        donor, timepoint, condition, replicate = parse_column_name(col_name)
        col_info[i] = (donor, timepoint, condition, replicate)
        if donor not in data_dict:
            data_dict[donor] = {}
        if timepoint not in data_dict[donor]:
            data_dict[donor][timepoint] = []
        data_dict[donor][timepoint].append(i)

# Read all data rows
with open(input_file, 'r', newline='') as infile:
    reader = csv.reader(infile, delimiter=delim)
    header = next(reader)  # skip re-read header
    data_lines = list(reader)

# Define condition ordering
condition_order = {'Mock': 0, 'FLU': 1}

# Sort columns for each donor/timepoint
for donor in data_dict:
    for timepoint in data_dict[donor]:
        # Extract columns for this donor/timepoint
        cols = data_dict[donor][timepoint]

        # Convert col indices to sortable tuples
        # We have col_info[i] = (donor, timepoint, condition, replicate)
        sortable_cols = []
        for c in cols:
            _, _, cond, rep = col_info[c]
            sortable_cols.append((condition_order.get(cond, 99), rep, c))

        # Sort by condition first, then by replicate number
        sortable_cols.sort(key=lambda x: (x[0], x[1]))

        # Update data_dict with sorted order of indices
        data_dict[donor][timepoint] = [x[2] for x in sortable_cols]

# Create directories and write out files
for donor in data_dict:
    donor_dir = os.path.join(output_dir, donor)
    if not os.path.exists(donor_dir):
        os.makedirs(donor_dir)

    for timepoint in data_dict[donor]:
        timepoint_dir = os.path.join(donor_dir, timepoint)
        if not os.path.exists(timepoint_dir):
            os.makedirs(timepoint_dir)

        output_csv = os.path.join(timepoint_dir, 'counts_named.csv')

        # Columns to include: Geneid + relevant columns (already sorted)
        col_indices = data_dict[donor][timepoint]
        all_indices = [geneid_index] + col_indices
        selected_header = [header[i] for i in all_indices]

        with open(output_csv, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter=delim)
            # Write header
            writer.writerow(selected_header)
            # Write rows
            for row in data_lines:
                selected_row = [row[i] for i in all_indices]
                writer.writerow(selected_row)
