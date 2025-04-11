import pandas as pd
import os
delim = ','

donor = 'res'
file_path= f'/home/pospim/Desktop/work/bioinformatics/proteomics/jitka/PXD040901/MQ_MERS/results/'
file = "dge_result.csv"
data = pd.read_csv(f'{file_path}{file}', delimiter=delim, engine='python')

print(data['Label'].head())

data['Label'] = data['Label'].astype(str)

mock_vs_variant = data[data['Label'].str.contains('Mock', na=False)]
mock_vs_variant_4h = mock_vs_variant[mock_vs_variant['Label'].str.contains('Mock_4h', na=False)]
mock_vs_variant_12h = mock_vs_variant[mock_vs_variant['Label'].str.contains('Mock_12h', na=False)]
mock_vs_variant_24h = mock_vs_variant[mock_vs_variant['Label'].str.contains('Mock_24h', na=False)]

output_dir = f'{file_path}res'
os.makedirs(output_dir, exist_ok=True)

variants_4h = mock_vs_variant_4h['Label'].str.extract(r'(\w+)_4h-Mock_4h')[0].dropna().unique()
variants_12h = mock_vs_variant_12h['Label'].str.extract(r'(\w+)_12h-Mock_12h')[0].dropna().unique()
variants_24h = mock_vs_variant_24h['Label'].str.extract(r'(\w+)_24h-Mock_24h')[0].dropna().unique()

variants = set(variants_4h) | set(variants_12h) | set(variants_24h)

for variant in variants:
    variant_dir = os.path.join(output_dir, variant)
    os.makedirs(variant_dir, exist_ok=True)

    # Filter data for the specific variant
    variant_data_4h = mock_vs_variant_4h[mock_vs_variant_4h['Label'].str.contains(f'{variant}_4h-Mock_4h', na=False)]
    variant_data_12h = mock_vs_variant_12h[mock_vs_variant_12h['Label'].str.contains(f'{variant}_12h-Mock_12h', na=False)]
    variant_data_24h = mock_vs_variant_24h[mock_vs_variant_24h['Label'].str.contains(f'{variant}_24h-Mock_24h', na=False)]


    # Save as CSV if data exists
    if not variant_data_4h.empty:
        print(f"Saving {variant}_4h.csv")
        variant_data_4h.to_csv(os.path.join(variant_dir, f'{variant}_4h.csv'), index=False, sep=delim)
    if not variant_data_12h.empty:
        print(f"Saving {variant}_12h.csv")
        variant_data_12h.to_csv(os.path.join(variant_dir, f'{variant}_12h.csv'), index=False, sep=delim)
    if not variant_data_24h.empty:
        print(f"Saving {variant}_24h.csv")
        variant_data_24h.to_csv(os.path.join(variant_dir, f'{variant}_24h.csv'), index=False, sep=delim)

print(f"Data successfully processed! Output saved in: {output_dir}")
