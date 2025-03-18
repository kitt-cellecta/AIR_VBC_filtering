import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from sklearn.mixture import GaussianMixture

def find_gmm_threshold(data):
    """Finds a threshold using Gaussian Mixture Model (GMM)."""
    if len(data) < 2:
        return None

    data_array = np.log10(data[data > 0].values).reshape(-1, 1)

    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data_array)

    mean1, mean2 = gmm.means_.flatten()
    var1, var2 = gmm.covariances_.flatten()

    if mean1 > mean2:
        mean1, mean2 = mean2, mean1
        var1, var2 = var2, var1

    a = 1 / (2 * var1) - 1 / (2 * var2)
    b = mean2 / var2 - mean1 / var1
    c = mean1**2 / (2 * var1) - mean2**2 / (2 * var2) + np.log(var2 / var1)
    intersection = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)

    return 10**intersection

def main(input_file, output_prefix):
    df = pd.read_csv(input_file, sep='\t')

    clone_barcode_counts = df.groupby('cloneId')['tagValueUMIBC'].nunique()

    clone_total_reads = df.groupby('cloneId')['readCount'].sum().reset_index()
    clone_total_reads['barcode_count'] = clone_total_reads['cloneId'].map(clone_barcode_counts)

    filtered_data = clone_total_reads[clone_total_reads['barcode_count'].between(1, 8)]

    thresholds = {}
    for barcode_count in range(1, 9):
        subset = filtered_data[filtered_data['barcode_count'] == barcode_count]['readCount']
        if not subset.empty:
            thresholds[barcode_count] = find_gmm_threshold(subset)

    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 16))
    axes = axes.flatten()

    for i, barcode_count in enumerate(range(1, 9)):
        subset = filtered_data[filtered_data['barcode_count'] == barcode_count]

        if not subset.empty:
            ax = axes[i]
            sns.histplot(subset['readCount'], bins=30, kde=True, ax=ax, log_scale=True)

            threshold = thresholds.get(barcode_count)
            if threshold:
                ax.axvline(threshold, color='red', linestyle='dashed', linewidth=2)
                ax.text(threshold, ax.get_ylim()[1] * 0.8, f'Threshold: {int(threshold)}',
                        color='red', ha='right', fontsize=10, weight='bold')

            ax.set_title(f'Clones in {barcode_count} Barcodes')
            ax.set_xlabel('Total Reads (log scale)')
            ax.set_ylabel('Frequency')

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_histograms.pdf")
    plt.close()

    filtered_data = filtered_data[
        filtered_data.apply(lambda row: row['readCount'] >= thresholds.get(row['barcode_count'], np.inf), axis=1)
    ]

    final_data = df[df['cloneId'].isin(filtered_data['cloneId'])].copy()

    grouped_final_data = final_data.groupby('cloneId').first().reset_index()
    grouped_final_data['readCount'] = final_data.groupby('cloneId')['readCount'].sum().values
    grouped_final_data['barcode_count'] = grouped_final_data['cloneId'].map(clone_barcode_counts)
    grouped_final_data = grouped_final_data.drop(columns=['tagValueUMIBC'], errors='ignore')

    output_file = f"{output_prefix}_filtered.tsv"
    grouped_final_data.sort_values("readCount", ascending=False).to_csv(output_file, sep='\t', index=False)
    print(f"Filtered data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter clones based on Gaussian Mixture Model thresholding.")
    parser.add_argument("input_file", type=str, help="Path to input TSV file.")
    parser.add_argument("output_prefix", type=str, help="Prefix for output files.")
    args = parser.parse_args()
    main(args.input_file, args.output_prefix)
