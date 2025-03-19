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
    parser.add_argument('percentage_threshold', help='Threshold value for filtering barcode hopping')
    args = parser.parse_args()
	
    run_filter_dir = os.path.dirname(os.path.abspath(__file__))
	
    if not os.path.isdir(args.dir):
        sys.exit(f"Error: Directory '{args.dir}' does not exist")

    try:
        report_file_name = f"{args.dir}VBC_Hopping_Filtering.report.txt"
        
        with open(report_file_name, 'w') as report_file:
            
            # Concatenate files
            output_all = os.path.join(args.dir, f"{args.sample}.clones_ALL.tsv")
            py1_concatenate = os.path.join(run_filter_dir, "concatenate_clones_files.py")
            run_subprocess(
                ["python3", py1_concatenate, args.dir],
                report_file
            )
            
            if not os.path.exists(output_all):
                sys.exit(f"Error: Failed to create {output_all}")
		
            # Filter Barcode Hopping 
            output_bcHop_filtered = os.path.join(args.dir, f"{args.sample}.clones_ALL.bcHop_filtered.tsv")
            py2_hopping = os.path.join(run_filter_dir, "barcode_hopping_filter.py")
            run_subprocess(
                ["python3", py2_hopping, output_all, args.percentage_threshold],
                report_file
            )
            
            if not os.path.exists(output_bcHop_filtered):
                sys.exit(f"Error: Failed to create {output_bcHop_filtered}")
		
            # Filter VBC
            output_vbc_filtered = os.path.join(args.dir, f"{args.sample}.clones_ALL.vbc_filtered.tsv")
            output_vbc_filtered_prefix = output_vbc_filtered.rsplit('.', 1)[0]
            py3_vbcFilt = os.path.join(run_filter_dir, "filter_vbc_2.py")
            run_subprocess(
                ["python3", py3_vbcFilt, output_bcHop_filtered, output_vbc_filtered_prefix],
                report_file
            )
            
            if not os.path.exists(output_vbc_filtered):
                sys.exit(f"Error: Stage 2 failed to create {output_vbc_filtered}")

            # Split filtered clones
            py4_split = os.path.join(run_filter_dir, "split_filtered_clones.py")
            run_subprocess(
                ["python3", py4_split, output_vbc_filtered],
                report_file
            )

            print("Barcode hopping and mutated sequence filtering completed...")
            os.remove(output_all)
            os.remove(output_bcHop_filtered)
            os.remove(output_vbc_filtered)
		
    except subprocess.CalledProcessError as e:
        sys.exit(f"Pipeline failed at {e.cmd[1]} with error {e.returncode}")

if __name__ == "__main__":
    main()