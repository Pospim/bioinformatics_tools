import pandas as pd
import os

separator = "\t"
# Set working directory
workdir = "/home/pospim/Desktop/work/bioinformatics/datasets/SARS/GSE213759_SARS_/Alpha/"

# Load CSV files
file_a = "counts_named_mock.csv"
file_b = "counts_named.csv"
df_a = pd.read_csv(os.path.join(workdir, file_a), sep=separator)
df_b = pd.read_csv(os.path.join(workdir, file_b), sep=separator)

print(df_a.head())
print(df_b.head())

# Merge dataframes on Geneid
df_merged = pd.merge(df_a, df_b, on="Geneid", how="inner")
df_merged.to_csv(os.path.join(workdir, "merged_counts.csv"), index=False, sep=separator)