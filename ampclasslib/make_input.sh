#!/usr/bin/env bash

# AmpliconClassifier input index builder
# Usage: make_input.sh <search_path1> [search_path2 ...] <output_name>
# arg 1-N: locations to search for AA files
# arg N+1: output name prefix

set -euo pipefail

# Validate arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <search_path1> [search_path2 ...] <output_name>" >&2
    exit 1
fi

# Extract output prefix (last argument) and search paths (all but last)
outpre="${@: -1}"
search_paths=("${@:1:$(($#-1))}")

# Temporary files
cycles_file="scf.txt"
graph_file="sgf.txt"
summary_file="ssf.txt"
sample_names="san.txt"
summary_names="ssn.txt"

# Clean up any existing temporary files
cleanup_temp_files() {
    rm -f "$cycles_file" "$graph_file" "$summary_file" "$sample_names" "$summary_names"
}

# Set up cleanup on exit
trap cleanup_temp_files EXIT

# Initialize temporary files
: > "$cycles_file"
: > "$graph_file"

# Find all AA graph and cycles files from AA runs
echo "Searching for AA cycles and graph files..."
for search_path in "${search_paths[@]}"; do
    exppath=$(realpath "$search_path")

    # Find cycles files with exclusions
    find "$exppath" -name "*_cycles.txt" \
        | grep -v "annotated_cycles" \
        | grep -v "/\._" \
        | grep -v "BPG_converted" \
        | grep -v "_classification/files/" \
        | sort >> "$cycles_file"

    # Find graph files with exclusions
    find "$exppath" -name "*_graph.txt" \
        | grep -v "features_to_graph" \
        | grep -v "/\._" \
        | grep -v "feature_to_graph" \
        | grep -v "_classification/files/" \
        | sort >> "$graph_file"
done

# Verify equal numbers of cycles and graph files
cycles_count=$(wc -l < "$cycles_file")
graph_count=$(wc -l < "$graph_file")

if [ "$cycles_count" -ne "$graph_count" ]; then
    echo "ERROR: Unequal numbers of cycles and graph files found! AA may not have completed correctly." >&2
    echo "Found $cycles_count cycles files and $graph_count graph files." >&2
    exit 1
fi

echo "Found $cycles_count matching cycles/graph file pairs"

# Extract sample names and create main input file
rev "$cycles_file" \
    | cut -f 1 -d '/' \
    | cut -c12- \
    | rev \
    | sed 's/_amplicon[0-9]*$//' > "$sample_names"

paste "$sample_names" "$cycles_file" "$graph_file" > "${outpre}.input"

# Find summary files for samples that may not have corresponding AA outputs
echo "Searching for summary files..."
: > "$summary_file"

for search_path in "${search_paths[@]}"; do
    exppath=$(realpath "$search_path")

    find "$exppath" -name "*_summary.txt" \
        | grep -v "/\._" \
        | grep -v "_classification/files/" \
        | sort >> "$summary_file"
done

summary_count=$(wc -l < "$summary_file")
echo "Found $summary_count summary files"

# Extract sample names from summary files and create summary map
rev "$summary_file" \
    | cut -f 1 -d '/' \
    | cut -c13- \
    | rev \
    | sed 's/_summary\.txt$//' > "$summary_names"

paste "$summary_names" "$summary_file" > "${outpre}_summary_map.txt"

echo "Created ${outpre}.input with $cycles_count entries"
echo "Created ${outpre}_summary_map.txt with $summary_count entries"

