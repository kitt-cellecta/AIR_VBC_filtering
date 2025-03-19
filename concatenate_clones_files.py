import os
import glob
import pandas as pd
import sys

print("Running concatenation...")

# Get the directory path from command line argument
if len(sys.argv) < 2:
    print("Please provide the directory path as an argument.")
    sys.exit(1)

directory = sys.argv[1]

# Find all TSV files matching the pattern
file_pattern = os.path.join(directory, '*clones_*.tsv')
files = glob.glob(file_pattern)

if not files:
    print("No matching files found.")
    sys.exit(1)

# Define TCR and BCR chains
tcr_chains = ['TRA', 'TRB', 'TRG', 'TRD']
bcr_chains = ['IGH', 'IGK', 'IGL']
all_chains = tcr_chains + bcr_chains

# Read all files and collect their column names
all_columns = set()
tcr_dataframes = []
bcr_dataframes = []
all_dataframes = []

sample = os.path.basename(files[0]).split('.clones_')[0]

for file in files:
    df = pd.read_csv(file, sep='\t', low_memory=False)
    
    # Extract the chain name from the filename
    filename = os.path.basename(file)
    chain = filename.split('clones_')[-1].split('.tsv')[0]
    
    # Add the 'chain' column
    df['chain'] = chain
    
    all_columns.update(df.columns)
    
    #if chain in tcr_chains:
    #    tcr_dataframes.append(df)
    #elif chain in bcr_chains:
    #    bcr_dataframes.append(df)
        
    all_dataframes.append(df)

# Ensure 'chain' is in all_columns
all_columns.add('chain')

# Function to process and save dataframes
def process_and_save(dataframes, output_file):

    if os.path.exists(output_file):
        print(f"File {output_file} already exists...")
        return False
    
    if not dataframes:
        print(f"No data for {output_file}")
        return
    
    # Reindex all dataframes with the complete set of columns
    for i, df in enumerate(dataframes):
        dataframes[i] = df.reindex(columns=list(all_columns), fill_value='')
    
    # Concatenate all dataframes
    combined_df = pd.concat(dataframes, ignore_index=True)
    
    # Save the combined dataframe to a new TSV file
    combined_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Combined file saved as {output_file}")
    print(f"Total number of clonotypes: {len(combined_df)}")

# Process and save TCR data
#tcr_output = os.path.join(directory, f'{sample}.clones_TCR.tsv')
#process_and_save(tcr_dataframes, tcr_output)

# Process and save BCR data
#bcr_output = os.path.join(directory, f'{sample}.clones_BCR.tsv')
#process_and_save(bcr_dataframes, bcr_output)

# Process and save ALL chains data
all_output = os.path.join(directory, f'{sample}.clones_ALL.tsv')
process_and_save(all_dataframes, all_output)

print()