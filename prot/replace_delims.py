import os
import pandas as pd

def replace_char_for_char(parent_dir, old_char, new_char):

    for root, dirs, files in os.walk(parent_dir):
        for file in files:
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)
                print(f"Processing {file_path}")

                try:
                    df = pd.read_csv(file_path, sep='\t', engine='python')
                    df.replace(old_char, new_char,regex=True, inplace=True)
                    df = df[~df['Protein'].str.startswith('SARS')]
                    df.to_csv(file_path, index=False, sep='\t')
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")


if __name__ == '__main__':
    parent_dir = '/home/pospim/Desktop/work/bioinformatics/proteomics/jitka/PXD036968/tmp'
    old_char = ';'
    new_char = '|'
    replace_char_for_char(parent_dir, old_char, new_char)