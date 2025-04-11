import pandas as pd

separator = '\t'
treshold = 8

workdir = "/home/pospim/Desktop/Work/bioinformatics/datasets/SARS/tst/"
merged_df = pd.read_csv(workdir + 'merged_deseq.csv', sep=separator)

def find_common_genes(df, treshold):
    log2_columns = [col for col in df.columns if 'Log2Fold' in col]
    df['gene_count'] = df[log2_columns].notna().sum(axis=1)

    common_genes_df = df[df['gene_count'] >= treshold]
    common_genes_df.to_csv(workdir + f'common_genes>={treshold}.csv', index=False, sep=separator)

def find_trends(df, treshold):
   # Identify all the log2FoldChange columns
    log2_columns = [col for col in df.columns if 'Log2FoldChange' in col]

    # Initialize empty columns to store results
    df['upregulated_count'] = (df[log2_columns] > 0).sum(axis=1)
    df['downregulated_count'] = (df[log2_columns] < 0).sum(axis=1)

    # Determine if the gene is predominantly upregulated or downregulated
    df['regulation_trend'] = df.apply(
        lambda row: 'Upregulated' if row['upregulated_count'] > row['downregulated_count']
        else 'Downregulated', axis=1
    )

    # Save the output to a new CSV file
    df.to_csv(workdir + f'gene_trend>={treshold}.csv', index=False, sep=separator)

    print("Analysis complete. Results saved to 'gene_regulation_trend.csv'.")