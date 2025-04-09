import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from sklearn.mixture import GaussianMixture
#from scipy.stats import norm
#from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks, argrelextrema
from sklearn.neighbors import KernelDensity

def find_gmm_threshold(data):
    """Finds a threshold using Gaussian Mixture Model (GMM)."""
    if len(data) < 2:
        return None

    data_array = np.log10(data[data > 0].values).reshape(-1, 1)

    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data_array)

    mean1, mean2 = gmm.means_.flatten()
    var1, var2 = gmm.covariances_.flatten()
    weight1, weight2 = gmm.weights_.flatten()

    if mean1 > mean2:
        mean1, mean2 = mean2, mean1
        var1, var2 = var2, var1
        weight1, weight2 = weight2, weight1

    a = 1 / (2 * var1) - 1 / (2 * var2)
    b = mean2 / var2 - mean1 / var1
    c = mean1**2 / (2 * var1) - mean2**2 / (2 * var2) + np.log(var2 / var1)
    intersection = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)

    return 10**intersection
    
    # Alternative solution: Minimum point detection method
    #def gaussian_mixture(x):
    #    return -1 * (weight1 * norm.pdf(x, mean1, np.sqrt(var1)) +
    #                 weight2 * norm.pdf(x, mean2, np.sqrt(var2)))
    #                 
    #result = minimize_scalar(gaussian_mixture, bracket=(mean1, mean2))
    #threshold = result.x
    #
    #return 10**threshold
    
def find_kde_mimima_threshold(data):
    """Finds threshold using KDE minima detection."""
    if len(data) < 2:
        return None
    
    data_array = np.log10(data[data > 0].values).reshape(-1, 1)
    
    # Adaptive bandwidth selection
    bandwidth = 0.5 * data_array.std()  # Silverman's rule of thumb factor
    if bandwidth == 0 or np.isnan(bandwidth):
        return None
    
    # Perform KDE analysis
    kde = KernelDensity(bandwidth=bandwidth)
    kde.fit(data_array)
    x = np.linspace(data_array.min(), data_array.max(), 1000).reshape(-1, 1)
    log_dens = kde.score_samples(x)
    
    # find locations of minima
    minima = argrelextrema(log_dens, np.less)[0]
    if len(minima) == 0:
        return None
    
    # Find first major minimum (ignore small fluctuations)
    main_minima = minima[log_dens[minima] < np.percentile(log_dens, 25)]
    if len(main_minima) == 0:
        return 10**x[minima[0]][0]  # Fallback to first minimum
    
    return 10**x[main_minima[0]][0]

def find_kde_mimima_threshold_2(data, barcode_count, output_prefix):
    """Finds threshold using KDE minima detection between two highest maxima."""
    if len(data) < 2:
        return None
    
    data_array = np.log10(data[data > 0].values).reshape(-1, 1)
    
    # Adaptive bandwidth selection
    bandwidth = 0.5 * data_array.std()
    if bandwidth == 0 or np.isnan(bandwidth):
        return None
    
    # Perform KDE analysis
    kde = KernelDensity(bandwidth=bandwidth)
    kde.fit(data_array)
    x = np.linspace(data_array.min(), data_array.max(), 1000).reshape(-1, 1)
    log_dens = kde.score_samples(x)
    
    # Find local minima and maxima
    minima = argrelextrema(log_dens, np.less)[0]
    maxima = argrelextrema(log_dens, np.greater)[0]
    
    # Check if we have at least two maxima
    if len(maxima) == 0:
        return 0
    
    if len(maxima) == 1:
        # Save maxima locations to a file
        with open(f"{output_prefix}.kde.maximas.txt", "a") as f:
            f.write(f"{barcode_count}\t{10**x[left_max][0]}\t{10**x[right_max][0]}\n")
		
		# Return fixed value as the cutoff if there's only one peak
        return 2
            
    if len(maxima) >= 2:
        # Get two highest maxima (by log density)
        sorted_maxima = maxima[np.argsort(-log_dens[maxima])]
        top_two_max = sorted_maxima[:2]
        left_max, right_max = sorted(top_two_max)
        
        # Save maxima locations to a file
        with open(f"{output_prefix}.kde.maximas.txt", "a") as f:
            f.write(f"{barcode_count}\t{10**x[left_max][0]}\t{10**x[right_max][0]}\n")
        
        # Find minima between these two maxima
        between_minima = minima[(minima > left_max) & (minima < right_max)]
        if between_minima.size > 0:
            # Select the most prominent minimum (lowest log density)
            selected_min = between_minima[np.argmin(log_dens[between_minima])]
            return 10**x[selected_min][0]
    
def main(input_file, output_prefix):
    
    print("Running VBC filtering...")
    
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    output_file = f"{output_prefix}.tsv"
    
    if os.path.exists(output_file):
        print(f"File {output_file} already exists...")
        return False
    
    # Calculate barcode counts and total reads per clone
    clone_barcode_counts = df.groupby('cloneId')['tagValueMIBC'].nunique()
    clone_total_reads = df.groupby('cloneId')['readCount'].sum().reset_index()
    clone_total_reads['barcode_count'] = clone_total_reads['cloneId'].map(clone_barcode_counts)

    # Perform KDE analysis
    thresholds = {}
    for barcode_count in range(1, 9):
        subset = clone_total_reads[clone_total_reads['barcode_count'] == barcode_count]['readCount']
        if not subset.empty:
            thresholds[barcode_count] = find_kde_mimima_threshold_2(subset, barcode_count, output_prefix)

    # Plot histograms and KDEs for clones
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))
    axes = axes.flatten()

    for i, barcode_count in enumerate(range(1, 9)):
        subset = clone_total_reads[clone_total_reads['barcode_count'] == barcode_count]

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

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.histograms.pdf")
    plt.close()

    # Apply thresholds to clones
    clone_total_reads = clone_total_reads.copy()
    clone_total_reads.loc[:, 'keep'] = clone_total_reads.apply(
        lambda row: row['readCount'] >= thresholds.get(row['barcode_count'], np.inf),
        axis=1
    )
    clones_to_keep = clone_total_reads[clone_total_reads['keep']]['cloneId']
    final_data = df[df['cloneId'].isin(clones_to_keep)].copy()

    # Group final data by cloneId and calculate summary statistics
    grouped_final_data = final_data.groupby('cloneId').first().reset_index()
    grouped_final_data['readCount'] = final_data.groupby('cloneId')['readCount'].sum().values
    grouped_final_data['barcode_count'] = grouped_final_data['cloneId'].map(clone_barcode_counts)
    grouped_final_data = grouped_final_data.drop(columns=['tagValueMIBC'], errors='ignore')

    # Save the final filtered data to a file
    grouped_final_data.sort_values("readCount", ascending=False).to_csv(output_file, sep='\t', index=False)
    print(f"Filtered data saved to {output_file}")

    # Print summary statistics
    original_rows = df.shape[0]
    filtered_rows = int(grouped_final_data['barcode_count'].sum())
    original_rows_readSum = int(df['readCount'].sum())
    filtered_rows_readSum = int(grouped_final_data['readCount'].sum())
    print(f"Original number of rows: {original_rows:,}")
    print(f"Original number of reads: {original_rows_readSum:,}")
    print(f"Final number of rows after mutated sequence filtering: {filtered_rows:,}")
    print(f"Final number of reads after mutated sequence filtering: {filtered_rows_readSum:,}")
    print(f"Reads removed: {original_rows_readSum - filtered_rows_readSum:,} "
          f"({(original_rows_readSum - filtered_rows_readSum)/original_rows_readSum:.1%})")
    print()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter clones based on Gaussian Mixture Model thresholding.")
    parser.add_argument("input_file", type=str, help="Path to input TSV file.")
    parser.add_argument("output_prefix", type=str, help="Prefix for output files.")
    args = parser.parse_args()
    main(args.input_file, args.output_prefix)
