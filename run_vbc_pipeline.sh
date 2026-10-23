#!/usr/bin/env bash
# Run VBC filtering followed by template-molecule quantification for one sample.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash run_vbc_pipeline.sh DIRECTORY SAMPLE_NAME [--mode MODE] [--peak-selection SELECTION]

Arguments:
  DIRECTORY                 MiXCR sample directory containing the clone table(s).
  SAMPLE_NAME               MiXCR output prefix for the sample.
  --mode MODE               bulk, single_cell, or bead (default: bulk).
  --peak-selection SELECTION
                            by_distance, by_prominence, or manual
                            (default: by_distance).

Set PYTHON_BIN to use a specific Python interpreter; otherwise python3 is used.
EOF
}

if [[ $# -lt 2 ]]; then
    usage >&2
    exit 2
fi

directory=$1
sample_name=$2
shift 2
mode=bulk
peak_selection=by_distance

while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)
            [[ $# -ge 2 ]] || { echo "Error: --mode requires a value." >&2; exit 2; }
            mode=$2
            shift 2
            ;;
        --peak-selection)
            [[ $# -ge 2 ]] || { echo "Error: --peak-selection requires a value." >&2; exit 2; }
            peak_selection=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unrecognized argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case $mode in bulk|single_cell|bead) ;; *) echo "Error: invalid mode: $mode" >&2; exit 2 ;; esac
case $peak_selection in by_distance|by_prominence|manual) ;; *) echo "Error: invalid peak selection: $peak_selection" >&2; exit 2 ;; esac

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON_BIN:-python3}

"$python_bin" "$script_dir/1_vbc_filtering/run_filter_vbc.py" \
    "$directory" "$sample_name" \
    --mode "$mode" \
    --peak-selection "$peak_selection"

if [[ $directory == */ ]]; then
    sample_dir=$directory
else
    sample_dir="$directory/"
fi
filtered_table="${sample_dir}${sample_name}.clones_ALL.filtered.tsv"
quant_report="${sample_dir}${sample_name}.templateEstimation.report.txt"
quantify_args=("$directory" "$sample_name" --mode "$mode")

# A manual rerun can replace the filtered table after prior quantification.
# Force a refresh only when that newly written input is newer than the report.
if [[ -f $filtered_table && -f $quant_report && $filtered_table -nt $quant_report ]]; then
    quantify_args+=(--force)
fi

"$python_bin" "$script_dir/2_quantify/run_quantify.py" "${quantify_args[@]}"
