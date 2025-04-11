import os, glob
import pandas as pd

delim = ','
def filter_csv_by_pval_log2FC(root_dir: str, pvalue_threshold: float=0.05, log2FC_threshold: float=2)->None:
    """
    Recursively searches for CSV files in `root_dir`, removes rows missing values
    in 'pvalue' or 'log2fold', filters rows where 'pvalue' < pval_threshold,
    and overwrites each CSV with the cleaned/filtered data.

    Parameters:
    -----------
    root_dir : str
        The top-level directory containing subdirectories with CSV files.
    pval_threshold : float
        The cutoff value for filtering rows based on 'pvalue'.
    """
    csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

    for csv_file in csv_files:
        df = pd.read_csv(csv_file, sep=delim, engine='python')
        required_columns = ['pvalue', 'log2FC']
        missing_columns = [col for col in required_columns if col not in df.columns]

        if missing_columns:
            print(f"Missing columns {missing_columns} in {csv_file}")
            continue
        df = df.dropna(subset=['pvalue', 'log2FC'])
        df = df[df['pvalue'] < pvalue_threshold]
        df = df[df['log2FC'].abs() > log2FC_threshold]
        if 'adj.pvalue' in df.columns:
            df.rename(columns={'adj.pvalue': 'padj'}, inplace=True)

        df = df.sort_values(by='pvalue', ascending=True)
        #df.rename(columns={'log2FC': 'log2fold'}, inplace=True)
        df.to_csv(csv_file, index=False, sep='\t')
        print(f"Filtered {csv_file}. Rows kept {len(df)}")

if __name__ == '__main__':
    root_dir = '/home/pospim/Desktop/work/bioinformatics/proteomics/jitka/PXD040901/MQ_MERS/results/MERS'
    filter_csv_by_pval_log2FC(root_dir)