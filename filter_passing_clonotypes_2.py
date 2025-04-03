import os
import sys
import pandas as pd

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
    df_other.to_csv("~/Desktop/tests/test.tsv", sep="\t", index=True)
	
    # Merge pivot table with metadata
    df_final = df_pivot.reset_index().merge(
        df_other,
        on='cloneId',
        how='left'
    )
    
    return df_final

def filter_passing_clonotypes(directory):
    if not os.path.isdir(directory):
        print(f"Error: The directory '{directory}' does not exist.")
        sys.exit(1)

    files = os.listdir(directory)
    
    # Find filtered clonotypes files
    file1 = next((f for f in files if f.endswith(".clones_ALL.bcHop_filtered.tsv")), None)
    file2 = next((f for f in files if f.endswith(".clones_ALL.vbc_filtered.tsv")), None)

    if not file1 or not file2:
        print("Error: Required files not found in the directory.")
        sys.exit(1)

    sample_name = file1.split(".clones_")[0]

    file2_path = os.path.join(directory, file2)
    df_file2 = pd.read_csv(file2_path, sep="\t")

    file1_path = os.path.join(directory, file1)
    df_file1 = pd.read_csv(file1_path, sep="\t")

    # Filter rows in barcode hopping filtered file where targetSequences match those in VBC filtered file
    unique_target_sequences = set(df_file2["targetSequences"].dropna())
    df_filtered = df_file1[df_file1["targetSequences"].isin(unique_target_sequences)]
    
	# Create simplified data frame with simplify_clonotypes_table
    df_simplified = simplify_clonotypes_table(df_filtered)

    # Group the filtered data by chain and save to separate files
    for chain, group in df_simplified.groupby("chain"):
        output_filename = f"{sample_name}.clones_{chain.upper()}_filtered.tsv"
        output_path = os.path.join(directory, output_filename)
        group = group.sort_values(by=['readCount'], ascending = False)
        group.to_csv(output_path, sep="\t", index=False)
        print(f"Filtered rows for chain {chain} saved to: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 filter_passing_clonotype_vbcs.py $DIRECTORY")
        sys.exit(1)

    input_directory = sys.argv[1]
    filter_passing_clonotypes(input_directory)