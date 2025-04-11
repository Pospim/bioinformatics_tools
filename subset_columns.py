import pandas as pd
import os

# Function to check if all substrings are in the column name
def contains_all_substrings(column_name, substrings):
    return all(substring in column_name for substring in substrings)

variant = ''
sample = ''
timepoint = '12h'
substrings = [sample, timepoint]

resdir = f'/home/pospim/Desktop/work/bioinformatics/datasets/SARS/GSE151513_/{variant}/{sample}/{timepoint}/'
workdir = f'/home/pospim/Desktop/work/bioinformatics/datasets/SARS/GSE151513_/{variant}/'

# Load dataset
df = pd.read_csv(workdir + 'counts_named.csv', encoding='UTF-8', sep='\t')
geneid_col = df[['Geneid']]


# Extract column containing all substrings
filtered_columns = [col for col in df.columns if contains_all_substrings(col, substrings)]
sub_df = df[filtered_columns]
print(sub_df)

# Perform left join
result_df = geneid_col.join(sub_df)

print(result_df)

# Save to new file
if result_df.shape[1] > 1:
    os.makedirs(resdir, exist_ok=True)
    result_df.to_csv(f'{resdir}counts_named.csv', index=False, encoding='UTF-8', sep='\t')
