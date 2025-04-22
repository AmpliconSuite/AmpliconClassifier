#!/bin/bash

# Default output filename
OUTPUT_FILE="cnv_bed_map.tsv"

print_help() {
    echo "Usage: $0 [--output FILE] [--help]"
    echo ""
    echo "Generates a TSV file mapping sample names to their absolute CNV BED file paths."
    echo ""
    echo "Options:"
    echo "  --output FILE   Specify the name of the output file (default: cnv_bed_map.tsv)"
    echo "  --help          Show this help message and exit"
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --help)
            print_help
            exit 0
            ;;
        --output)
            shift
            OUTPUT_FILE="$1"
            ;;
        *)
            echo "Unknown option: $1"
            print_help
            exit 1
            ;;
    esac
    shift
done

echo "Creating CNV BED file mapping..."
echo "Output will be written to: $OUTPUT_FILE"

# Find and process *_CNV_CALLS.bed files
cnv_files=$(find . -type f -name "*_CNV_CALLS.bed")
if [ -z "$cnv_files" ]; then
    echo "No CNV BED files found matching '*_CNV_CALLS.bed'. Exiting."
    exit 1
fi

# Overwrite the output file
> "$OUTPUT_FILE"

count=0
while IFS= read -r filepath; do
    filename=$(basename "$filepath")
    sample_name=$(echo "$filename" | sed -E 's/(.*)_CNV_CALLS\.bed/\1/' | sed 's/\.cs\.rmdup//')
    abs_path=$(realpath "$filepath")
    echo -e "${sample_name}\t${abs_path}" >> "$OUTPUT_FILE"
    ((count++))
done <<< "$cnv_files"

echo "Done! Found $count CNV BED files."
echo "Mapping saved to: $OUTPUT_FILE"
