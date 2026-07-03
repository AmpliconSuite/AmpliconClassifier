#!/usr/bin/env bash

# AmpliconClassifier input index builder
# Usage: make_input.sh <search_path1> [search_path2 ...] <output_name>
# arg 1-N: locations to search for AA files
# arg N+1: output name prefix

set -u

usage() {
    echo "Usage: $0 <search_path1> [search_path2 ...] <output_name>" >&2
}

if [ "$#" -lt 2 ]; then
    usage
    exit 1
fi

outpre="${@: -1}"
search_paths=()
for ((i=1; i<$#; i++)); do
    search_paths+=("${!i}")
done

outdir=$(dirname "$outpre")
if [ "$outdir" != "." ] && [ ! -d "$outdir" ]; then
    mkdir -p "$outdir" || {
        echo "ERROR: Could not create output directory: $outdir" >&2
        exit 1
    }
fi

resolved_paths=()
for search_path in "${search_paths[@]}"; do
    if ! exppath=$(realpath "$search_path" 2>/dev/null); then
        echo "Warning: Cannot resolve path $search_path, skipping..." >&2
        continue
    fi

    if [ ! -d "$exppath" ]; then
        echo "Warning: Search path is not a directory: $exppath, skipping..." >&2
        continue
    fi

    resolved_paths+=("$exppath")
done

if [ "${#resolved_paths[@]}" -eq 0 ]; then
    echo "ERROR: No valid search paths provided." >&2
    exit 1
fi

is_excluded_path() {
    case "$1" in
        */._*|*BPG_converted*|*_classification/files/*)
            return 0
            ;;
    esac

    return 1
}

sample_name_from_cycles() {
    local base
    base=$(basename "$1")
    base=${base%_cycles.txt}
    echo "$base" | sed 's/_amplicon[0-9]*$//'
}

sample_name_from_summary() {
    local base
    base=$(basename "$1")
    base=${base%_summary.txt}
    # CoRAL writes "{sample}_amplicon_summary.txt" while AA writes "{sample}_summary.txt".
    # Strip a trailing "_amplicon" so the summary sample name matches the cycles-derived name
    # (sample_name_from_cycles strips "_amplicon[0-9]*"). No-op for AA-style names.
    base=${base%_amplicon}
    echo "$base"
}

cycles_files=()
graph_files=()
summary_files=()

echo "Searching for AA cycles and graph files..."
for exppath in "${resolved_paths[@]}"; do
    while IFS= read -r file; do
        if is_excluded_path "$file"; then
            continue
        fi
        if [[ "$(basename "$file")" == *annotated_cycles* ]]; then
            continue
        fi
        cycles_files+=("$file")
    done < <(find "$exppath" -type f -name "*_cycles.txt" 2>/dev/null | sort)

    while IFS= read -r file; do
        if is_excluded_path "$file"; then
            continue
        fi
        case "$(basename "$file")" in
            *features_to_graph*|*feature_to_graph*)
                continue
                ;;
        esac
        graph_files+=("$file")
    done < <(find "$exppath" -type f -name "*_graph.txt" 2>/dev/null | sort)

    while IFS= read -r file; do
        if is_excluded_path "$file"; then
            continue
        fi
        summary_files+=("$file")
    done < <(find "$exppath" -type f -name "*_summary.txt" 2>/dev/null | sort)
done

declare -A graph_seen
for graph_file in "${graph_files[@]}"; do
    graph_seen["$graph_file"]=0
done

missing_graph_count=0
: > "${outpre}.input" || {
    echo "ERROR: Could not write ${outpre}.input" >&2
    exit 1
}

for cycles_file in "${cycles_files[@]}"; do
    graph_file="${cycles_file%_cycles.txt}_graph.txt"
    if [ ! -f "$graph_file" ]; then
        echo "ERROR: Missing graph file for cycles file:" >&2
        echo "  cycles: $cycles_file" >&2
        echo "  expected graph: $graph_file" >&2
        missing_graph_count=$((missing_graph_count + 1))
        continue
    fi

    graph_seen["$graph_file"]=1
    sample_name=$(sample_name_from_cycles "$cycles_file")
    printf "%s\t%s\t%s\n" "$sample_name" "$cycles_file" "$graph_file" >> "${outpre}.input"
done

orphan_graph_count=0
for graph_file in "${graph_files[@]}"; do
    if [ "${graph_seen[$graph_file]}" -eq 0 ]; then
        echo "Warning: Graph file has no matching cycles file, ignoring: $graph_file" >&2
        orphan_graph_count=$((orphan_graph_count + 1))
    fi
done

if [ "$missing_graph_count" -gt 0 ]; then
    echo "ERROR: Found $missing_graph_count cycles file(s) without matching graph files." >&2
    exit 1
fi

input_count=$(wc -l < "${outpre}.input" 2>/dev/null || echo 0)
if [ "$input_count" -eq 0 ]; then
    echo "No cycles/graph file pairs found - created empty input file"
else
    echo "Found $input_count matching cycles/graph file pairs"
fi

if [ "$orphan_graph_count" -gt 0 ]; then
    echo "Ignored $orphan_graph_count graph file(s) without matching cycles files" >&2
fi

echo "Searching for summary files..."
: > "${outpre}_summary_map.txt" || {
    echo "ERROR: Could not write ${outpre}_summary_map.txt" >&2
    exit 1
}

for summary_file in "${summary_files[@]}"; do
    sample_name=$(sample_name_from_summary "$summary_file")
    printf "%s\t%s\n" "$sample_name" "$summary_file" >> "${outpre}_summary_map.txt"
done

summary_count=$(wc -l < "${outpre}_summary_map.txt" 2>/dev/null || echo 0)
if [ "$summary_count" -eq 0 ]; then
    echo "No summary files found - created empty summary map"
else
    echo "Found $summary_count summary files"
fi

echo "Created ${outpre}.input with $input_count entries"
echo "Created ${outpre}_summary_map.txt with $summary_count entries"
