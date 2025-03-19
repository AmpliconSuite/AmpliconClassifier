#!/bin/bash

# Set the output file name
OUTPUT_FILE="cnv_bed_map.tsv"

# Remove the output file if it already exists
if [ -f "$OUTPUT_FILE" ]; then
    rm "$OUTPUT_FILE"
fi

# Find all files matching the pattern and process them
find . -type f -name "*_CNV_CALLS.bed" | while read -r filepath; do
    # Get just the filename without the path
    filename=$(basename "$filepath")
    
    # Extract the part of the filename before the _CNV_*CALLS.bed pattern
    sample_name=$(echo "$filename" | sed -E 's/(.*)_CNV_CALLS\.bed/\1/')
    
    # Get the absolute path
    abs_path=$(realpath "$filepath")
    
    # Write to the output file
    echo -e "${sample_name}\t${abs_path}" >> "$OUTPUT_FILE"
done

echo "Processing complete. Results saved to $OUTPUT_FILE with $(wc -l < "$OUTPUT_FILE") entries"
