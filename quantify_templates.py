import pandas as pd
import os
import math
import argparse

def quantify_templates(directory, sample):
    
    # Check if the analysis has already been completed
    report_file = f"{sample}.templateEstimation.report.txt"
    report_path = os.path.join(directory, report_file)
    if os.path.exists(report_path):
        print(f"Template estimation for {sample} has already completed. Will not run again...")
        return
    
    # Get normFactor from kde.maximas file
    maxima_file = f"{sample}.clones_ALL.vbc_filtered.kde.maximas.txt"
    maxima_path = os.path.join(directory, maxima_file)
    
    try:
        maxima_df = pd.read_csv(maxima_path, sep='\t', header=None)
        normFactor = float(maxima_df.iloc[0, 2])
    except Exception as e:
        print(f"Error reading maxima file: {str(e)}")
        return

    # Process each chain type
    chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']
    
    processed_files = []
    for chain in chains:
        input_file = f"{sample}.clones_{chain}_filtered.tsv"
        input_path = os.path.join(directory, input_file)
        
        if not os.path.exists(input_path):
            continue
            
        try:

            df = pd.read_csv(input_path, sep='\t')            
            if 'readCount' in df.columns:
                df['templateEstimate'] = df['readCount'].apply(
                    lambda x: math.ceil(x / normFactor)
                )
                
                output_file = f"{sample}.clones_{chain}_quantified.tsv"
                output_path = os.path.join(directory, output_file)
                df.to_csv(output_path, sep='\t', index=False)
                processed_files.append(output_file)
                
            else:
                print(f"'readCount' column missing in {input_file}")
                
        except Exception as e:
            print(f"Error processing {chain}: {str(e)}")
            
    # Create report file
    with open(report_path, 'w') as report:
        report.write("Template Estimation Analysis Completed\n")
        report.write(f"Norm Factor: {normFactor}\n")
        report.write("Processed Files:\n")
        for file in processed_files:
            report.write(f"- {file}\n")
    
    print(f"Template estimation analysis completed for {sample}\n")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Quantify templates based on normFactor")
    parser.add_argument("directory", help="Directory containing the sample files")
    parser.add_argument("sample", help="Sample name")

    # Parse arguments
    args = parser.parse_args()
    
    # Run the processing function
    quantify_templates(args.directory, args.sample)

