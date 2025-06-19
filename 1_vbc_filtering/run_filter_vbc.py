#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import io

def run_subprocess(command, output_file):
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True
    )

    # Read and write output line by line
    for line in process.stdout:
        sys.stdout.write(line)  # Write to console
        output_file.write(line) # Write to file
        sys.stdout.flush()
        output_file.flush()

    # Wait for the process to complete
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)

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
        report_file_name = f"{directory}{sample_name}.vbcFiltering.report.txt"
        if os.path.exists(report_file_name):
            print(f"Clones have already been filtered for {sample_name}. See '{report_file_name}'...")
            sys.exit(0)
        
        with open(report_file_name, 'w') as report_file:
            
            # Concatenate files
            output_all = os.path.join(directory, f"{sample_name}.clones_ALL.tsv")
            py1_concatenate = os.path.join(run_filter_dir, "concatenate_clones_files.py")
            run_subprocess(
                ["python3", py1_concatenate, directory],
                report_file
            )
            
            if not os.path.exists(output_all):
                sys.exit(f"Error: Failed to create {output_all}")
		
            # QC Filtering steps with VBC
            output_filtered = os.path.join(directory, f"{sample_name}.clones_ALL.filtered.tsv")
            py2_filtering = os.path.join(run_filter_dir, "filter.py")
            run_subprocess(
                ["python3", py2_filtering, output_all, directory, sample_name],
                report_file
            )
            
            if not os.path.exists(output_filtered):
                sys.exit(f"Error: Failed to create {output_filtered}")
            
            print("Barcode hopping and mutated sequence filtering completed...")
            os.remove(output_all)
		
    except subprocess.CalledProcessError as e:
        sys.exit(f"Pipeline failed at {e.cmd[1]} with error {e.returncode}")

if __name__ == "__main__":
    main()