#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def filter_clones(input_file, percentage_threshold):
    
    print("Running barcode hopping filtering...")    
    
    # Parse arguments
    parser = argparse.ArgumentParser(description='Filter clone data by read count thresholds')
    parser.add_argument('input_file', help='Path to input TSV file (format: *.clones_ALL.tsv)')
    parser.add_argument('percentage', type=float, help='Percentage threshold for filtering (e.g., 5 for 5%)')
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    if args.percentage <= 0 or args.percentage >= 100:
        print("Error: Percentage must be between 0 and 100", file=sys.stderr)
        sys.exit(1)
    
    # Read TSV file
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    
    # Generate output filename
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}.bcHop_filtered.tsv"
    
    if os.path.exists(output_file):
        print(f"File {output_file} already exists...")
        return False
    
    # Calculate maximum read count per clone group
    df['group_max'] = df.groupby('cloneId')['readCount'].transform('max')
    
    # Filter rows where readCount >= fraction_threshold of group maximum
    fraction_threshold = percentage_threshold/100
    filtered_df = df[df['readCount'] >= fraction_threshold * df['group_max']].copy()
    
    # Clean up temporary column
    filtered_df.drop(columns=['group_max'], inplace=True)
    
    # Save filtered data
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"Filtered data saved to: {output_file}")
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter clone data by read count thresholds')
    parser.add_argument('input_file', help='Path to input TSV file (format: *.clones_ALL.tsv)')
    parser.add_argument('percentage', type=float, help='Percentage threshold for filtering (e.g., 5 for 5%)')
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    if args.percentage <= 0 or args.percentage >= 100:
        print("Error: Percentage must be between 0 and 100", file=sys.stderr)
        sys.exit(1)
    
    filter_clones(args.input_file, args.percentage)