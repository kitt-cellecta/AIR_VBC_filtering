#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os

def main():
    
    parser = argparse.ArgumentParser(description='Process clones files and filter low abundance clones')
    parser.add_argument('dir', help='Input directory containing initial data files')
    parser.add_argument('sample', help='Sample name prefix for output files')
    parser.add_argument('percentage_threshold', help='Threshold value for filtering barcode hopping')
    args = parser.parse_args()

    if not os.path.isdir(args.dir):
        sys.exit(f"Error: Directory '{args.dir}' does not exist")

    try:
        
        # Concatenate files
        output_all = os.path.join(args.dir, f"{args.sample}.clones_ALL.tsv")
        print("Running concatenation...")
        subprocess.run(
            ["python3", "concatenate_clones_files.py", args.dir],
            check=True
        )
        print()
        
        if not os.path.exists(output_all):
            sys.exit(f"Error: Failed to create {output_all}")
		
		# Filter Barcode Hopping 
        output_bcHop_filtered = os.path.join(args.dir, f"{args.sample}.clones_ALL.bcHop_filtered.tsv")
        print("Running barcode hopping filtering...")
        subprocess.run(
            ["python3", "barcode_hopping_filter.py", output_all, args.percentage_threshold],
            check=True
        )
        
        if not os.path.exists(output_all):
            sys.exit(f"Error: Failed to create {output_bcHop_filtered}")
        print()
		
        # Filter VBC
        output_vbc_filtered = os.path.join(args.dir, f"{args.sample}.clones_ALL.vbc_filtered.tsv")
        output_vbc_filtered_prefix = output_vbc_filtered.rsplit('.', 1)[0]
        print("Running VBC filtering...")
        subprocess.run(
            ["python3", "filter_vbc_2.py", output_bcHop_filtered, output_vbc_filtered_prefix],
            check=True
        )
        print()
        
        if not os.path.exists(output_vbc_filtered):
            sys.exit(f"Error: Stage 2 failed to create {output_vbc_filtered}")


        # Split filtered clones
        print("Splitting filtered clones...")
        subprocess.run(
            ["python3", "split_filtered_clones.py", output_vbc_filtered],
            check=True
        )
        print()

        print("Barcode hopping and mutated sequence filtering completed...")
        os.remove(output_all)
        os.remove(output_bcHop_filtered)
        os.remove(output_vbc_filtered)
		
    except subprocess.CalledProcessError as e:
        sys.exit(f"Pipeline failed at {e.cmd[1]} with error {e.returncode}")

if __name__ == "__main__":
    main()
