import argparse
import pandas as pd
import os
import sys

DELIM='\t'

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Map GSM/SRR column names in counts file to sample names using a mapping file."
    )
    parser.add_argument(
        "--counts",
        required=True,
        help="Path to the counts.csv file. Must contain a gene_id column and GSM/SRR column names.",
    )
    parser.add_argument(
        "--samples",
        required=True,
        help="Path to the samples.csv file. Must contain columns: sample, gsm, srr.",
    )
    parser.add_argument(
        "--id_column",
        default="gene_id",
        help="Name of the identifier column (default: gene_id).",
    )
    parser.add_argument(
        "--output",
        default="counts_named.csv",
        help="Output path for the mapped counts CSV file.",
    )
    return parser.parse_args()

def load_data(counts_file, samples_file, id_column):
    if not os.path.exists(counts_file):
        sys.exit(f"Counts file not found: {counts_file}")

    if not os.path.exists(samples_file):
        sys.exit(f"Samples file not found: {samples_file}")

    counts_df = pd.read_csv(counts_file, sep=DELIM)

    print(counts_df.head)

    if id_column not in counts_df.columns:
        sys.exit(f"ID column '{id_column}' not found in counts file.")

    samples_df = pd.read_csv(samples_file, sep=DELIM)

    required_cols = {"sample", "GSM", "SRR"}

    if not required_cols.issubset(samples_df.columns):
        sys.exit(f"Samples file must contain columns: {required_cols}")

    return counts_df, samples_df

def map_columns(counts_df, samples_df, id_column):
    # Extract GSM and SRR mappings
    gsm_map = dict(zip(samples_df["GSM"], samples_df["sample"]))
    srr_map = dict(zip(samples_df["SRR"], samples_df["sample"]))

    # Column mapping logic (excluding id_column)
    new_columns = []
    for col in counts_df.columns:
        if col == id_column:
            new_columns.append(col)
        elif col in gsm_map:
            new_columns.append(gsm_map[col])
        elif col in srr_map:
            new_columns.append(srr_map[col])
        else:
            new_columns.append(col)

    # Update columns in counts_df
    counts_df.columns = new_columns
    return counts_df

def save_mapped_counts(counts_df, output_file):
    counts_df.to_csv(output_file, index=False, sep=DELIM)
    print(f"✅ Mapped counts file saved to: {output_file}")

def main():
    args = parse_arguments()

    # Load files
    counts_df, samples_df = load_data(args.counts, args.samples, args.id_column)

    # Map columns
    mapped_counts_df = map_columns(counts_df, samples_df, args.id_column)

    # Save output
    save_mapped_counts(mapped_counts_df, args.output)

if __name__ == "__main__":
    main()