import os
import sys
import pandas as pd

def filter_passing_clonotypes(directory):
    if not os.path.isdir(directory):
        print(f"Error: The directory '{directory}' does not exist.")
        sys.exit(1)

    files = os.listdir(directory)
    
    # Find the required files
    file1 = next((f for f in files if f.endswith(".clones_ALL.bcHop_filtered.tsv")), None)
    file2 = next((f for f in files if f.endswith(".clones_ALL.vbc_filtered.tsv")), None)

    if not file1 or not file2:
        print("Error: Required files not found in the directory.")
        sys.exit(1)

    sample_name = file1.split(".clones_")[0]

    # Load data from file2
    file2_path = os.path.join(directory, file2)
    df_file2 = pd.read_csv(file2_path, sep="\t")

    # Extract unique targetSequences from file2
    unique_target_sequences = set(df_file2["targetSequences"].dropna())

    # Load data from file1
    file1_path = os.path.join(directory, file1)
    df_file1 = pd.read_csv(file1_path, sep="\t")

    # Filter rows in file1 where targetSequences match those in file2
    df_filtered = df_file1[df_file1["targetSequences"].isin(unique_target_sequences)]

    # Group the filtered data by chain and save to separate files
    for chain, group in df_filtered.groupby("chain"):
        output_filename = f"{sample_name}.clones_{chain.upper()}_filtered.tsv"
        output_path = os.path.join(directory, output_filename)
        group.to_csv(output_path, sep="\t", index=False)
        print(f"Filtered rows for chain {chain} saved to: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 filter_passing_clonotype_vbcs.py $DIRECTORY")
        sys.exit(1)

    input_directory = sys.argv[1]
    filter_passing_clonotypes(input_directory)