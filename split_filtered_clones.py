import os
import pandas as pd
import sys

def split_tsv_by_chain(input_file):
    if not os.path.isfile(input_file):
        print(f"Error: {input_file} is not a valid file.")
        return

    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # Get the sample name and directory
    directory = os.path.dirname(input_file)
    sample_name = os.path.basename(input_file).split('.clones_')[0]
    
    # Group by 'chain' column and save to separate files
    for chain, group in df.groupby('chain'):
        output_file = os.path.join(directory, f"{sample_name}.clones_{chain}_filtered.tsv")
        group.to_csv(output_file, sep='\t', index=False)
        print(f"Created file: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py /path/to/input_file.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    split_tsv_by_chain(input_file)