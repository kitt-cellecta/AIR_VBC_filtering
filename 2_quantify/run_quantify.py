#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import io
import pandas as pd

######
### main 
### inputs:
###     sample_name =  sample name
###     directory = location of mixcr files (set as "", i.e. blank, if working in Platforma)
### outputs: 
###     $SAMPLE_NAME.clones_$CHAIN.quantified.tsv - clonotype tables with molecular estimates
######

def main():
    
    parser = argparse.ArgumentParser(description='Process clones files and filter low abundance clones')
    parser.add_argument('dir', help='Input directory containing initial data files')
    parser.add_argument('sample', help='Sample name prefix for output files')
    args = parser.parse_args()
    sample_name = args.sample
    directory = args.dir
	
    if not directory.endswith(os.path.sep):
        directory += os.path.sep
	
    run_filter_dir = os.path.dirname(os.path.abspath(__file__))
	
    if not os.path.isdir(directory):
        sys.exit(f"Error: Directory '{directory}' does not exist")

    try:
            
        # run quantification process
        
        quantification_output = os.path.join(directory, f"{sample_name}.clones_ALL.quantified.tsv")
        py1_quantify = os.path.join(run_filter_dir, "normalize.py")
        subprocess.run(["python3", py1_quantify, directory, sample_name])
            
        if not os.path.exists(quantification_output):
            sys.exit(f"Error: Failed to create {quantification_output}")
		
		# split the quantified files into each unique chain value 
		
        df = pd.read_csv(quantification_output, sep='\t')
        chains = df['chain'].unique()
        for chain_value in chains:
            chain_df = df[df['chain'] == chain_value]
            output_file = os.path.join(directory, f"{sample_name}.clones_{chain_value}.quantified.tsv")
            chain_df.to_csv(output_file, sep='\t', index=False)
		
		# check if any of the chain specific files are generated
		
        chain_values = ['TRA', 'TRB', 'TRG', 'TRD', 'IGH', 'IGK', 'IGL']
        any_exists = any(
            os.path.exists(os.path.join(directory, f"{sample_name}.clones_{chain}.quantified.tsv"))
            for chain in chain_values
        )
        if any_exists and os.path.exists(quantification_output):
            os.remove(quantification_output)
            print(f"Removed input file: {quantification_output}")
        else:
            sys.exit("Cannot generate any of the chain-specific quantified clonotype tables...")
		
    except subprocess.CalledProcessError as e:
        sys.exit(f"Pipeline failed at {e.cmd[1]} with error {e.returncode}")

if __name__ == "__main__":
    main()