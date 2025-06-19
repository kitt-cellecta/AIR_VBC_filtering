#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity
from scipy.signal import find_peaks, argrelextrema
from scipy.stats import iqr

######
### barcode_hopping_filter -- barcode hopping filter; remove VBCs with low reads relative to other VBCs of the same clonotype sequence
###	inputs: 
###    input_file = file containing all clonotypes from MiXCR with VBC preset
###    percentage = cutoff percentage of reads that a VBC must have to be considered non-barcode hopping 
### outputs:
###    output = barcode hopping filtered tsv
######

def barcode_hopping_filter(input_file, percentage):
    
    print("Running barcode hopping filtering...")    
        
    # Validate Inputs
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    
    if percentage <= 0 or percentage >= 100:
        print("Error: Percentage must be between 0 and 100", file=sys.stderr)
        sys.exit(1)
    
    # Read TSV file
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    
    # Calculate maximum read count per clone group
    df['group_max'] = df.groupby('cloneId')['readCount'].transform('max')
    
    # Filter rows where readCount >= fraction_threshold of group maximum
    fraction_threshold = percentage/100
    filtered_df = df[df['readCount'] >= fraction_threshold * df['group_max']].copy()
    
    # Clean up temporary column
    filtered_df.drop(columns=['group_max'], inplace=True)
    
    # Print Stats
    original_rows = df.shape[0]
    filtered_rows = filtered_df.shape[0]    
    original_rows_readSum = int(df['readCount'].sum())
    filtered_rows_readSum = int(filtered_df['readCount'].sum())
    print(f"Original number of rows: {original_rows:,}")
    print(f"Original number of reads: {original_rows_readSum:,}")
    print(f"Final number of rows after barcode hopping filtering: {filtered_rows:,}")
    print(f"Final number of reads after barcode hopping filtering: {filtered_rows_readSum:,}")
    print(f"Reads removed: {original_rows_readSum - filtered_rows_readSum:,} "
          f"({(original_rows_readSum - filtered_rows_readSum)/original_rows_readSum:.1%})")
    print()

    # Return filtered data    
    return filtered_df

######
### find_kde_mimima_threshold -- function to model histograms by KDE and find maximas/minimas;
###    helper function for reads_per_clonotype_filter
######

def find_kde_mimima_threshold(data, barcode_count, sample_name):
    """Finds threshold using KDE minima detection between two highest maxima."""
    if len(data) < 2:
        return None
    
    data_array = np.log10(data[data > 0].values).reshape(-1, 1)
    
    # Adaptive bandwidth selection
    #bandwidth = 0.5 * data_array.std()
    n = len(data_array)
    sigma = np.std(data_array, ddof=1)
    bandwidth = 1.06 * sigma * (n ** (-1/5))
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
        return (2, 0, 0)
    
    if len(maxima) == 1:
        left_max = 0 # setting the first peak as non-existent at 0
        right_max = maxima[0]
		
        # Find minima between these two maxima
        between_minima = minima[(minima > left_max) & (minima < right_max)]
        if between_minima.size > 0:
            # Select the most prominent minimum (lowest log density)
            selected_min = between_minima[np.argmin(log_dens[between_minima])]
            return (10**x[selected_min][0], 10**x[left_max][0], 10**x[right_max][0])
            
    if len(maxima) >= 2:
        # Get two highest maxima (by log density)
        sorted_maxima = maxima[np.argsort(-log_dens[maxima])]
        top_two_max = sorted_maxima[:2]
        left_max, right_max = sorted(top_two_max)
                
        # Find minima between these two maxima
        between_minima = minima[(minima > left_max) & (minima < right_max)]
        if between_minima.size > 0:
            # Select the most prominent minimum (lowest log density)
            selected_min = between_minima[np.argmin(log_dens[between_minima])]
            return (10**x[selected_min][0], 10**x[left_max][0], 10**x[right_max][0])

######
### reads_per_clonotype_filter -- filter clonotypes with low number of reads using VBC bins
### inputs:
###     input_file = barcode hopping filter data frame
###     sample_name = prefix name of output file
### outputs: 
###     output_file = clonotypes filtered by the number of reads
######
    
def reads_per_clonotype_filter(input_file, directory, sample_name):
    
    print("Running VBC filtering...")
    
    df = input_file
    
    # Calculate barcode counts and total reads per clone
    clone_barcode_counts = df.groupby('cloneId')['tagValueMIBC'].nunique()
    clone_total_reads = df.groupby('cloneId')['readCount'].sum().reset_index()
    clone_total_reads['barcode_count'] = clone_total_reads['cloneId'].map(clone_barcode_counts)

    # Perform KDE analysis
    thresholds = {}
    left_maxes = {}
    right_maxes = {}
    for barcode_count in range(1, 9):
        subset = clone_total_reads[clone_total_reads['barcode_count'] == barcode_count]['readCount']
        if not subset.empty:
            thresholds[barcode_count], left_maxes[barcode_count], right_maxes[barcode_count] = find_kde_mimima_threshold(subset, barcode_count, sample_name)
    
    # Check results of KDE analysis across 8 barcodes before writing to file and performing filtering
    # in certain cases, the VBC = 4 or higher have no first peak but there might be a small peak after the second peak 
    # for a couple of clonotypes. we know that the increase of values from VBC = 1, 2, etc is a slight doubling, so
    # we're using this information to check if the value of the filter went haywire. the cutoff that we're using is 10x
    # increase.
    
    sorted_keys = sorted(thresholds.keys())
    
    for i in range(1, len(sorted_keys)):

        current_key = sorted_keys[i]
        previous_key = sorted_keys[i - 1]

        # check if previous threshold is 10x less than the current threshold (previousThreshold10x)
        previousThreshold10x = thresholds[current_key] > 10 * thresholds[previous_key]    

        # check if next threshold is 10x more than the current threshold (nextThreshold10x)   
        if current_key < len(sorted_keys):
            next_key = sorted_keys[i + 1]
            nextThreshold10x = thresholds[next_key] > 10 * thresholds[current_key]
        else:
            nextThreshold10x = False
        
        # check if the previous threshold is just a default threshold setting (previousThresholdDefault)
        previousThresholdDefault = thresholds[previous_key] == 2

        # modify thresholds file if previousThreshold10x, nextThreshold10x and previousThresholdDefault are all True
        if previousThreshold10x and nextThreshold10x and previousThresholdDefault: 
            right_maxes[current_key] = left_maxes[current_key] 	# the first peak is actually the second peak
            left_maxes[current_key] = 0							# the first peak is set at 0
            thresholds[current_key] = 2							# default cutoff value = 2
            
	# Save maxima and threshold locations to a file
    maximas_file = os.path.join(directory, f"{sample_name}.kde.maximas.txt")
    with open(maximas_file, "a") as f:
        for barcode_count in range(1, 9):
            f.write(f"{barcode_count}\t{thresholds[barcode_count]}\t{left_maxes[barcode_count]}\t{right_maxes[barcode_count]}\n")
    
    # Plot histograms and KDEs for clones
    histogram_file = os.path.join(directory, f"{sample_name}.histograms.pdf")
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
    plt.savefig(histogram_file)
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
    grouped_final_data = grouped_final_data.sort_values("readCount", ascending=False)

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
    
    return grouped_final_data

######
### simplify_clonotypes_table -- simplifies clonotype table by creating one row for each clonotype and extra columns for each VBC count;
###     helper function for filter_passing_clonotypes
######

def simplify_clonotypes_table(df_filtered):
    """
    Simplifies the clonotypes table by creating read counts column for each VBC 
    This reduces the number of rows on the table
    Then it creates a new column for the total read counts called readCount
    """
    # Create simplified table for readCount values, set value to 0 if the BC does not have a particular clonotype
    df_pivot = df_filtered.pivot_table(
        index='cloneId',
        columns='tagValueMIBC',
        values='readCount',
        fill_value=0
    ).add_prefix('readCount_')

    for bc in [f'BC{i}' for i in range(1, 9)]:
        col_name = f'readCount_{bc}'
        if col_name not in df_pivot.columns:
            df_pivot[col_name] = 0
    
    # Reorder columns and add total readCount column
    bc_columns = [f'readCount_BC{i}' for i in range(1, 9)]
    df_pivot = df_pivot[bc_columns]
    df_pivot['readCount'] = df_pivot.sum(axis=1)

    # Get metadata from row with max readCount per cloneId
    idx = df_filtered.groupby('cloneId')['readCount'].idxmax()
    other_columns = df_filtered.columns.difference(
        ['tagValueMIBC', 'readCount', 'readFraction']
    )
    df_other = df_filtered.loc[idx, other_columns].reset_index(drop=True)
	
    # Merge pivot table with metadata
    df_final = df_pivot.reset_index().merge(
        df_other,
        on='cloneId',
        how='left'
    )
    
    return df_final

######
### filter_passing_clonotypes -- generate final table of passing clonotypes
### inputs:
###     df_bcHop_filtered = barcode hopping filtered df (output of barcode_hopping_filter)
###     df_rpclon_filtered = reads per clonotype filtered df (output of reads_per_clonotype_filter)
###     sample_name = sample name
### outputs: 
###     df_passing_clonotypes = simplified df post-barcode hopping and reads per clonotype filters
######

def filter_passing_clonotypes(df_bcHop_filtered, df_rpclon_filtered):

    df_file1 = df_bcHop_filtered
    df_file2 = df_rpclon_filtered

    # Filter rows in barcode hopping filtered file where targetSequences match those in VBC filtered file
    unique_target_sequences = set(df_file2["targetSequences"].dropna())
    df_filtered = df_file1[df_file1["targetSequences"].isin(unique_target_sequences)]
    
	# Create simplified data frame with simplify_clonotypes_table
    df_simplified = simplify_clonotypes_table(df_filtered)

    return df_simplified
    
######
### main 
### inputs:
###     input_file = file containing all clonotypes from MiXCR with VBC preset
###     sample_name = sample name
### outputs: 
###     df_passing_clonotypes = simplified df post-barcode hopping and reads per clonotype filters
######
    
def main(input_file, directory, sample_name):
    """
    Step 1: barcode_hopping_filter -- filter barcode hopping 
    Step 2: reads_per_clonotype_filter -- filter low read count clonotypes
    Step 3: filter_passing_clonotypes -- generate final table
    """
	
    percentage = 5 # default setting
    mixcr_clones_file = input_file
	
    df_bcHop_filtered = barcode_hopping_filter(mixcr_clones_file, percentage)
    df_rpclon_filtered = reads_per_clonotype_filter(df_bcHop_filtered, directory, sample_name)
    df_passing_clonotypes = filter_passing_clonotypes(df_bcHop_filtered, df_rpclon_filtered)
    
    output_filename = f"{sample_name}.clones_ALL.filtered.tsv"
    output_path = os.path.join(directory, output_filename)
    df_passing_clonotypes = df_passing_clonotypes.sort_values(by=['readCount'], ascending = False)
    df_passing_clonotypes.to_csv(output_path, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quality control filtering using VBCs for Cellecta DriverMap AIR.")
    parser.add_argument("input_file", type=str, help="Path to input TSV file.")
    parser.add_argument("directory", type=str, help="Directory of MiXCR output files")
    parser.add_argument("sample_name", type=str, help="Prefix for output files.")
    args = parser.parse_args()
    main(args.input_file, args.directory, args.sample_name)
