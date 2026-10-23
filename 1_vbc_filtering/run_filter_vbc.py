#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import io
import csv
import json
import shutil


MANUAL_TEMPLATE_NAME = "peak_manual_select.csv"
MANIFEST_SUFFIX = ".vbc_analysis.manifest.json"


def manual_choices_are_all_keep(directory, mode):
    template = os.path.join(directory, MANUAL_TEMPLATE_NAME)
    if not os.path.isfile(template):
        raise FileNotFoundError(f"Manual selection template not found: {template}")
    with open(template, newline="") as handle:
        rows = list(csv.DictReader(handle))
    choices = {int(row["vbc_bin"]): row["filter_choice"].strip() for row in rows}
    if set(choices) != set(range(1, 9)):
        raise ValueError("Manual template must contain exactly VBC bins 1 through 8.")
    active_bins = range(1, 5) if mode == "single_cell" else range(1, 9)
    return all(choices[barcode_count] == "keep" for barcode_count in active_bins)


def archive_original_analysis(directory, sample_name, mode):
    archive_dir = os.path.join(directory, "original_vbc_analysis")
    # A prior manual attempt may have moved the automatic outputs before being
    # interrupted.  In that case, resume from the archived manifest rather
    # than requiring a manifest in the sample root.
    if os.path.exists(archive_dir):
        archived_manifest = os.path.join(archive_dir, f"{sample_name}{MANIFEST_SUFFIX}")
        if not os.path.isfile(archived_manifest):
            raise FileExistsError(f"Original VBC analysis already archived: {archive_dir}")
        with open(archived_manifest) as handle:
            metadata = json.load(handle)
        if metadata.get("status") != "completed" or metadata.get("manual"):
            raise ValueError("Manual mode requires a completed automatic VBC analysis.")
        if metadata.get("mode") != mode:
            raise ValueError(f"Manual mode expects original mode {metadata.get('mode')!r}, not {mode!r}.")
        return archive_dir

    manifest = os.path.join(directory, f"{sample_name}{MANIFEST_SUFFIX}")
    if not os.path.isfile(manifest):
        raise FileNotFoundError("Manual mode requires a completed original VBC analysis manifest.")
    with open(manifest) as handle:
        metadata = json.load(handle)
    if metadata.get("status") != "completed" or metadata.get("manual"):
        raise ValueError("Manual mode requires a completed automatic VBC analysis.")
    if metadata.get("mode") != mode:
        raise ValueError(f"Manual mode expects original mode {metadata.get('mode')!r}, not {mode!r}.")
    os.makedirs(archive_dir)
    names = [
        f"{sample_name}.clones_ALL.filtered.tsv",
        f"{sample_name}.histograms.pdf",
        f"{sample_name}.kde.maximas.txt",
        f"{sample_name}{MANIFEST_SUFFIX}",
        f"{sample_name}.vbcFiltering.report.txt",
    ]
    for name in names:
        source = os.path.join(directory, name)
        if os.path.exists(source):
            shutil.move(source, os.path.join(archive_dir, name))
    return archive_dir

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
    parser.add_argument('--mode', type=str, choices=['bulk', 'single_cell'], default='bulk', help="Processing mode: 'bulk' or 'single_cell' DriverMap AIR.")
    parser.add_argument(
        '--filter-selection', '--peak-selection',
        dest='filter_selection',
        choices=['by_distance', 'by_prominence', 'manual'],
        default='by_distance',
        help=(
            "Strategy for KDE distributions with more than two peaks "
            "(default: by_distance, the original behavior)."
        ),
    )
    args = parser.parse_args()
    sample_name = args.sample
    directory = args.dir
    mode = args.mode
    filter_selection = args.filter_selection
    if not directory.endswith(os.path.sep):
        directory += os.path.sep
    run_filter_dir = os.path.dirname(os.path.abspath(__file__))
    if not os.path.isdir(directory):
        sys.exit(f"Error: Directory '{directory}' does not exist")
    try:
        if filter_selection == "manual":
            if mode == "single_cell":
                print("Manual mode: VBC bins 5-8 are ignored in this 4-VBC mode.")
            if manual_choices_are_all_keep(directory, mode):
                print("Manual mode: all active VBC bins are 'keep'; no filtering changes requested.")
                return
            archive_original_analysis(directory, sample_name, mode)
        report_file_name = f"{directory}{sample_name}.vbcFiltering.report.txt"
        completed_manual_manifest = os.path.join(directory, f"{sample_name}{MANIFEST_SUFFIX}")
        if filter_selection == "manual" and os.path.isfile(completed_manual_manifest):
            with open(completed_manual_manifest) as handle:
                completed_metadata = json.load(handle)
            if completed_metadata.get("status") == "completed" and completed_metadata.get("manual"):
                print(f"Manual filtering has already completed for {sample_name}. See '{report_file_name}'...")
                return
        if filter_selection != "manual" and os.path.exists(report_file_name):
            print(f"Clones have already been filtered for {sample_name}. See '{report_file_name}'...")
            sys.exit(0)
        with open(report_file_name, 'w') as report_file:
            output_all = os.path.join(directory, f"{sample_name}.clones_ALL.tsv")
            output_all_preexisting = os.path.exists(output_all)
            py1_concatenate = os.path.join(run_filter_dir, "concatenate_clones_files.py")
            run_subprocess(
                [sys.executable, py1_concatenate, directory],
                report_file
            )
            if not os.path.exists(output_all):
                sys.exit(f"Error: Failed to create {output_all}")
            filter_input = output_all
            # QC Filtering steps with VBC
            output_filtered = os.path.join(directory, f"{sample_name}.clones_ALL.filtered.tsv")
            py2_filtering = os.path.join(run_filter_dir, "filter.py")
            filter_cmd = [sys.executable, py2_filtering, filter_input, directory, sample_name]
            if mode:
                filter_cmd += ["--mode", mode]
            filter_cmd += ["--filter-selection", filter_selection]
            run_subprocess(
                filter_cmd,
                report_file
            )
            if not os.path.exists(output_filtered):
                sys.exit(f"Error: Failed to create {output_filtered}")
            print("Barcode hopping and mutated sequence filtering completed...")
            if not output_all_preexisting:
                os.remove(output_all)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Pipeline failed at {e.cmd[1]} with error {e.returncode}")

if __name__ == "__main__":
    main()
