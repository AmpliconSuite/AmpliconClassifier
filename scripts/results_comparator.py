#!/usr/bin/env python3
"""
AmpliconClassifier Results Comparison Tool

Compares two AmpliconClassifier results tables to identify gains, losses,
classification changes, and structural rearrangements (splits/merges) of focal amplifications.
"""

import argparse
import pandas as pd
import sys
from collections import Counter, defaultdict
import re
from typing import List, Tuple, Dict, Set, Optional


def parse_coordinates(coord_str: str) -> List[Tuple[str, int, int]]:
    """
    Parse coordinate string like ['chr6:20350615-22839372', 'chr3:9521655-10184493']
    Returns list of (chromosome, start, end) tuples.
    """
    if pd.isna(coord_str) or coord_str == "['']" or coord_str == "" or coord_str == "[]":
        return []

    # Remove outer brackets and split by commas
    coord_str = coord_str.strip('[]')
    coords = []

    # Find all coordinate patterns - handle both quoted and unquoted
    pattern = r"['\"]?([^'\",:]+:[^'\",:]+)['\"]?"
    matches = re.findall(pattern, coord_str)

    for match in matches:
        if ':' in match and '-' in match:
            try:
                chrom, range_part = match.split(':')
                start, end = map(int, range_part.split('-'))
                coords.append((chrom, start, end))
            except ValueError:
                continue

    return coords


def calculate_overlap(coords1: List[Tuple[str, int, int]], coords2: List[Tuple[str, int, int]]) -> int:
    """
    Calculate total overlap in base pairs between two sets of coordinates.
    """
    total_overlap = 0

    for chrom1, start1, end1 in coords1:
        for chrom2, start2, end2 in coords2:
            if chrom1 == chrom2:
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                if overlap_end > overlap_start:
                    total_overlap += overlap_end - overlap_start

    return total_overlap


def is_valid_feature(feature: pd.Series, run: int) -> bool:
    """
    Check if a feature represents a real focal amplification.
    Returns False for samples that were analyzed but had no features.
    Prints warning for features with classification but empty coordinates.
    """
    # Check if coordinates are empty or invalid
    coords = parse_coordinates(str(feature.get('Location', '')))
    classification = str(feature.get('Classification', '')).strip()
    feature_id = str(feature.get('Feature ID', ''))

    if not coords:
        # Check if this is a real feature with empty coordinates (warning case)
        if classification and classification not in ['', 'nan', 'None'] and not feature_id.endswith('_NA'):
            print(
                f"WARNING: Input {run} feature {feature_id} has classification '{classification}' but empty coordinates. Ignoring this feature.")
        return False

    # Check for placeholder Feature IDs that indicate no features
    if '_NA' in feature_id or feature_id.endswith('_NA'):
        return False

    return True


def features_overlap(feature1: pd.Series, feature2: pd.Series, min_overlap: int = 1000) -> bool:
    """
    Check if two features overlap by at least min_overlap base pairs.
    """
    # First check if both features are valid
    if not is_valid_feature(feature1, 1) or not is_valid_feature(feature2, 2):
        return False

    coords1 = parse_coordinates(feature1['Location'])
    coords2 = parse_coordinates(feature2['Location'])

    if not coords1 or not coords2:
        return False

    overlap = calculate_overlap(coords1, coords2)
    return overlap >= min_overlap


def parse_gene_list(gene_str: str) -> List[str]:
    """
    Parse gene list string like ['CCND3', 'TFEB'] into list of genes.
    """
    if pd.isna(gene_str) or gene_str == "['']" or gene_str == "" or gene_str == "[]":
        return []

    # Remove outer brackets and split by comma
    gene_str = gene_str.strip('[]')

    # Find all gene names - handle both quoted and unquoted
    pattern = r"['\"]?([^'\",:]+)['\"]?"
    genes = re.findall(pattern, gene_str)

    # Filter out empty strings and clean whitespace
    return [gene.strip() for gene in genes if gene.strip()]


def get_representative_genes(oncogenes_lists: List[List[str]], all_genes_lists: List[List[str]], max_genes: int = 3) -> \
        List[str]:
    """
    Get representative genes that are common across overlapping features.
    Prioritize oncogenes, fall back to all genes if no oncogenes.
    """
    # Filter out empty lists and ensure all items are strings
    clean_oncogenes = []
    for genes in oncogenes_lists:
        if genes and isinstance(genes, list):
            clean_genes = [str(g) for g in genes if pd.notna(g) and str(g).strip()]
            if clean_genes:
                clean_oncogenes.append(clean_genes)

    clean_all_genes = []
    for genes in all_genes_lists:
        if genes and isinstance(genes, list):
            clean_genes = [str(g) for g in genes if pd.notna(g) and str(g).strip()]
            if clean_genes:
                clean_all_genes.append(clean_genes)

    # Find intersection of oncogenes
    if clean_oncogenes:
        oncogene_sets = [set(genes) for genes in clean_oncogenes]
        if oncogene_sets:
            common_oncogenes = set.intersection(*oncogene_sets)
            if common_oncogenes:
                return sorted(list(common_oncogenes))[:max_genes]

    # Find intersection of all genes
    if clean_all_genes:
        gene_sets = [set(genes) for genes in clean_all_genes]
        if gene_sets:
            common_genes = set.intersection(*gene_sets)
            if common_genes:
                return sorted(list(common_genes))[:max_genes]

    # If no common genes, take union and pick representative ones
    all_oncogenes = set()
    all_all_genes = set()

    for genes in clean_oncogenes:
        all_oncogenes.update(genes)
    for genes in clean_all_genes:
        all_all_genes.update(genes)

    if all_oncogenes:
        return sorted(list(all_oncogenes))[:max_genes]
    elif all_all_genes:
        return sorted(list(all_all_genes))[:max_genes]
    else:
        return ["No genes"]


def calculate_jaccard_index(coords1: List[Tuple[str, int, int]], coords2: List[Tuple[str, int, int]]) -> float:
    """
    Calculate Jaccard index between two coordinate sets at base-pair level.
    Returns intersection_bp / union_bp.
    """
    if not coords1 or not coords2:
        return 0.0

    # Convert coordinates to intervals per chromosome
    chrom_intervals1 = defaultdict(list)
    chrom_intervals2 = defaultdict(list)

    for chrom, start, end in coords1:
        chrom_intervals1[chrom].append((start, end))

    for chrom, start, end in coords2:
        chrom_intervals2[chrom].append((start, end))

    total_intersection = 0
    total_union = 0

    # Process each chromosome
    all_chroms = set(chrom_intervals1.keys()) | set(chrom_intervals2.keys())

    for chrom in all_chroms:
        intervals1 = chrom_intervals1.get(chrom, [])
        intervals2 = chrom_intervals2.get(chrom, [])

        # Merge overlapping intervals within each set
        merged1 = merge_intervals(intervals1)
        merged2 = merge_intervals(intervals2)

        # Calculate intersection and union for this chromosome
        intersection_bp = calculate_interval_intersection(merged1, merged2)
        union_bp = calculate_interval_union(merged1, merged2)

        total_intersection += intersection_bp
        total_union += union_bp

    return total_intersection / total_union if total_union > 0 else 0.0


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping intervals."""
    if not intervals:
        return []

    sorted_intervals = sorted(intervals)
    merged = [sorted_intervals[0]]

    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged


def calculate_interval_intersection(intervals1: List[Tuple[int, int]], intervals2: List[Tuple[int, int]]) -> int:
    """Calculate total base pairs of intersection between two interval sets."""
    total_intersection = 0

    for start1, end1 in intervals1:
        for start2, end2 in intervals2:
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            if overlap_end > overlap_start:
                total_intersection += overlap_end - overlap_start

    return total_intersection


def calculate_interval_union(intervals1: List[Tuple[int, int]], intervals2: List[Tuple[int, int]]) -> int:
    """Calculate total base pairs of union between two interval sets."""
    # Combine all intervals and merge
    all_intervals = intervals1 + intervals2
    merged = merge_intervals(all_intervals)

    return sum(end - start for start, end in merged)


def format_coordinates(coords: List[Tuple[str, int, int]]) -> str:
    """
    Format coordinates for display.
    """
    if not coords:
        return "Unknown"

    coord_strs = [f"{chrom}:{start}-{end}" for chrom, start, end in coords]
    return ", ".join(coord_strs)


def create_bipartite_mapping(df1: pd.DataFrame, df2: pd.DataFrame) -> Dict[str, Dict[str, List[str]]]:
    """
    Create bipartite mapping between features in two runs.
    Returns dict with sample names as keys, and mapping info as values.
    """
    mapping = {}

    # Get all samples present in either run
    samples1 = set(df1['Sample name'].unique())
    samples2 = set(df2['Sample name'].unique())
    all_samples = samples1.union(samples2)

    for sample in all_samples:
        sample_df1 = df1[df1['Sample name'] == sample] if sample in samples1 else pd.DataFrame()
        sample_df2 = df2[df2['Sample name'] == sample] if sample in samples2 else pd.DataFrame()

        # Filter out invalid features (no-feature samples)
        if not sample_df1.empty:
            sample_df1 = sample_df1[sample_df1.apply(lambda row: is_valid_feature(row, 1), axis=1)]
        if not sample_df2.empty:
            sample_df2 = sample_df2[sample_df2.apply(lambda row: is_valid_feature(row, 2), axis=1)]

        # Create overlap matrix
        overlaps = defaultdict(list)

        for idx1, feature1 in sample_df1.iterrows():
            for idx2, feature2 in sample_df2.iterrows():
                if features_overlap(feature1, feature2):
                    overlaps[f"run1_{idx1}"].append(f"run2_{idx2}")
                    overlaps[f"run2_{idx2}"].append(f"run1_{idx1}")

        # Add unmatched features
        for idx1, _ in sample_df1.iterrows():
            if f"run1_{idx1}" not in overlaps:
                overlaps[f"run1_{idx1}"] = []

        for idx2, _ in sample_df2.iterrows():
            if f"run2_{idx2}" not in overlaps:
                overlaps[f"run2_{idx2}"] = []

        mapping[sample] = {
            'overlaps': dict(overlaps),
            'df1': sample_df1,
            'df2': sample_df2
        }

    return mapping


def analyze_events(mapping: Dict) -> List[Dict]:
    """
    Analyze the bipartite mapping to identify different types of events.
    """
    events = []

    for sample, data in mapping.items():
        overlaps = data['overlaps']
        df1 = data['df1']
        df2 = data['df2']

        # Group features by their overlap patterns
        processed_features = set()

        # Process each connected component
        for feature_id in overlaps:
            if feature_id in processed_features:
                continue

            # Find all connected features (BFS)
            to_visit = [feature_id]
            component = set()

            while to_visit:
                current = to_visit.pop(0)
                if current in component:
                    continue
                component.add(current)

                for connected in overlaps[current]:
                    if connected not in component:
                        to_visit.append(connected)

            # Mark all features in this component as processed
            processed_features.update(component)

            # Analyze this component
            run1_features = [f for f in component if f.startswith('run1_')]
            run2_features = [f for f in component if f.startswith('run2_')]

            event = analyze_component(sample, run1_features, run2_features, df1, df2)
            if event:
                events.append(event)

    return events


def analyze_component(sample: str, run1_features: List[str], run2_features: List[str],
                      df1: pd.DataFrame, df2: pd.DataFrame) -> Optional[Dict]:
    """
    Analyze a connected component to determine the type of event.
    """
    run1_count = len(run1_features)
    run2_count = len(run2_features)

    # Get feature data
    run1_data = []
    run2_data = []

    for feat_id in run1_features:
        idx = int(feat_id.split('_')[1])
        if idx in df1.index:
            run1_data.append(df1.loc[idx])

    for feat_id in run2_features:
        idx = int(feat_id.split('_')[1])
        if idx in df2.index:
            run2_data.append(df2.loc[idx])

    if run1_count == 0 and run2_count == 0:
        return None

    # Determine event type
    if run1_count == 1 and run2_count == 0:
        event_type = "LOSS"
    elif run1_count == 0 and run2_count == 1:
        event_type = "GAIN"
    elif run1_count == 1 and run2_count == 1:
        # Check for classification change
        if run1_data[0]['Classification'] != run2_data[0]['Classification']:
            event_type = "CLASS_CHANGE"
        else:
            event_type = "STABLE"
    elif run1_count == 1 and run2_count > 1:
        event_type = "SPLIT"
    elif run1_count > 1 and run2_count == 1:
        event_type = "MERGE"
    else:
        event_type = "COMPLEX_MERGE_SPLIT"

    # Collect gene information
    run1_oncogenes = []
    run1_all_genes = []
    run2_oncogenes = []
    run2_all_genes = []

    for f in run1_data:
        oncogenes = parse_gene_list(str(f.get('Oncogenes', '')))
        all_genes = parse_gene_list(str(f.get('All genes', '')))
        run1_oncogenes.append(oncogenes)
        run1_all_genes.append(all_genes)

    for f in run2_data:
        oncogenes = parse_gene_list(str(f.get('Oncogenes', '')))
        all_genes = parse_gene_list(str(f.get('All genes', '')))
        run2_oncogenes.append(oncogenes)
        run2_all_genes.append(all_genes)

    all_oncogenes = run1_oncogenes + run2_oncogenes
    all_all_genes = run1_all_genes + run2_all_genes

    representative_genes = get_representative_genes(all_oncogenes, all_all_genes)

    # Get coordinates
    run1_coords = []
    run2_coords = []

    for f in run1_data:
        run1_coords.extend(parse_coordinates(f.get('Location', '')))
    for f in run2_data:
        run2_coords.extend(parse_coordinates(f.get('Location', '')))

    location = format_coordinates(run1_coords + run2_coords)

    # Calculate Jaccard overlap for applicable event types
    jaccard_overlap = None
    if event_type in ["STABLE", "CLASS_CHANGE", "SPLIT", "MERGE", "COMPLEX_MERGE_SPLIT"]:
        jaccard_overlap = calculate_jaccard_index(run1_coords, run2_coords)

    # Get copy numbers and classifications
    run1_cn = []
    run2_cn = []
    run1_classes = []
    run2_classes = []

    for f in run1_data:
        cn = f.get('Feature median copy number', 0)
        # Handle potential NaN or string values
        try:
            cn = float(cn) if pd.notna(cn) else 0.0
        except (ValueError, TypeError):
            cn = 0.0
        run1_cn.append(cn)
        run1_classes.append(str(f.get('Classification', '')))

    for f in run2_data:
        cn = f.get('Feature median copy number', 0)
        # Handle potential NaN or string values
        try:
            cn = float(cn) if pd.notna(cn) else 0.0
        except (ValueError, TypeError):
            cn = 0.0
        run2_cn.append(cn)
        run2_classes.append(str(f.get('Classification', '')))

    return {
        'Sample': sample,
        'Event_Type': event_type,
        'Run1_Features': [str(f.get('Feature ID', '')) for f in run1_data],
        'Run2_Features': [str(f.get('Feature ID', '')) for f in run2_data],
        'Run1_Classifications': run1_classes,
        'Run2_Classifications': run2_classes,
        'Run1_Copy_Numbers': run1_cn,
        'Run2_Copy_Numbers': run2_cn,
        'Location': location,
        'Representative_Genes': representative_genes,
        'Jaccard_Overlap': jaccard_overlap
    }


def write_hierarchical_report(events: List[Dict], output_file: str):
    """
    Write hierarchical summary report.
    """
    # Group events by sample and type
    sample_events = defaultdict(lambda: defaultdict(list))

    for event in events:
        sample = event['Sample']
        event_type = event['Event_Type']
        sample_events[sample][event_type].append(event)

    with open(output_file, 'w') as f:
        f.write("AmpliconClassifier Comparison Report\n")
        f.write("=" * 50 + "\n\n")

        if not sample_events:
            f.write("No differences found between the two runs.\n")
            return

        for sample in sorted(sample_events.keys()):
            f.write(f"SAMPLE: {sample}\n")
            f.write("-" * (len(sample) + 8) + "\n\n")

            events_for_sample = sample_events[sample]

            # Gains
            if 'GAIN' in events_for_sample:
                f.write(f"  Gains ({len(events_for_sample['GAIN'])}):\n")
                for event in events_for_sample['GAIN']:
                    classes = ", ".join(event['Run2_Classifications'])
                    cns = ", ".join([f"{cn:.1f}" for cn in event['Run2_Copy_Numbers']])
                    genes = ", ".join(event['Representative_Genes'])
                    f.write(f"    - {classes}: {event['Location']} (CN={cns}, {genes})\n")
                f.write("\n")

            # Losses
            if 'LOSS' in events_for_sample:
                f.write(f"  Losses ({len(events_for_sample['LOSS'])}):\n")
                for event in events_for_sample['LOSS']:
                    classes = ", ".join(event['Run1_Classifications'])
                    cns = ", ".join([f"{cn:.1f}" for cn in event['Run1_Copy_Numbers']])
                    genes = ", ".join(event['Representative_Genes'])
                    f.write(f"    - {classes}: {event['Location']} (CN={cns}, {genes})\n")
                f.write("\n")

            # Classification Changes
            if 'CLASS_CHANGE' in events_for_sample:
                f.write(f"  Classification Changes ({len(events_for_sample['CLASS_CHANGE'])}):\n")
                for event in events_for_sample['CLASS_CHANGE']:
                    old_class = event['Run1_Classifications'][0]
                    new_class = event['Run2_Classifications'][0]
                    old_cn = event['Run1_Copy_Numbers'][0]
                    new_cn = event['Run2_Copy_Numbers'][0]
                    genes = ", ".join(event['Representative_Genes'])
                    jaccard_str = f", Jaccard={event['Jaccard_Overlap']:.3f}" if event[
                                                                                     'Jaccard_Overlap'] is not None else ""
                    f.write(
                        f"    - {event['Location']}: {old_class}→{new_class} (CN: {old_cn:.1f}→{new_cn:.1f}, {genes}{jaccard_str})\n")
                f.write("\n")

            # Splits
            if 'SPLIT' in events_for_sample:
                f.write(f"  Splits ({len(events_for_sample['SPLIT'])}):\n")
                for event in events_for_sample['SPLIT']:
                    old_class = event['Run1_Classifications'][0]
                    old_cn = event['Run1_Copy_Numbers'][0]
                    genes = ", ".join(event['Representative_Genes'])
                    jaccard_str = f", Jaccard={event['Jaccard_Overlap']:.3f}" if event[
                                                                                     'Jaccard_Overlap'] is not None else ""
                    f.write(f"    - {old_class}: {event['Location']} (CN={old_cn:.1f}, {genes}{jaccard_str}) →\n")
                    for i, (new_class, new_cn) in enumerate(
                            zip(event['Run2_Classifications'], event['Run2_Copy_Numbers'])):
                        f.write(f"        {new_class} (CN={new_cn:.1f})\n")
                f.write("\n")

            # Merges
            if 'MERGE' in events_for_sample:
                f.write(f"  Merges ({len(events_for_sample['MERGE'])}):\n")
                for event in events_for_sample['MERGE']:
                    new_class = event['Run2_Classifications'][0]
                    new_cn = event['Run2_Copy_Numbers'][0]
                    genes = ", ".join(event['Representative_Genes'])
                    jaccard_str = f", Jaccard={event['Jaccard_Overlap']:.3f}" if event[
                                                                                     'Jaccard_Overlap'] is not None else ""
                    f.write(
                        f"    - Multiple features → {new_class}: {event['Location']} (CN={new_cn:.1f}, {genes}{jaccard_str})\n")
                    for old_class, old_cn in zip(event['Run1_Classifications'], event['Run1_Copy_Numbers']):
                        f.write(f"        ← {old_class} (CN={old_cn:.1f})\n")
                f.write("\n")

            # Stable events
            if 'STABLE' in events_for_sample:
                f.write(f"  Stable ({len(events_for_sample['STABLE'])}):\n")
                for event in events_for_sample['STABLE']:
                    classes = event['Run1_Classifications'][0]  # Same for both runs
                    old_cn = event['Run1_Copy_Numbers'][0]
                    new_cn = event['Run2_Copy_Numbers'][0]
                    genes = ", ".join(event['Representative_Genes'])
                    jaccard_str = f", Jaccard={event['Jaccard_Overlap']:.3f}" if event[
                                                                                     'Jaccard_Overlap'] is not None else ""
                    f.write(
                        f"    - {classes}: {event['Location']} (CN: {old_cn:.1f}→{new_cn:.1f}, {genes}{jaccard_str})\n")
                f.write("\n")

            # Complex merge-splits
            if 'COMPLEX_MERGE_SPLIT' in events_for_sample:
                f.write(f"  Complex Merge-Splits ({len(events_for_sample['COMPLEX_MERGE_SPLIT'])}):\n")
                for event in events_for_sample['COMPLEX_MERGE_SPLIT']:
                    genes = ", ".join(event['Representative_Genes'])
                    jaccard_str = f", Jaccard={event['Jaccard_Overlap']:.3f}" if event[
                                                                                     'Jaccard_Overlap'] is not None else ""
                    f.write(f"    - Complex rearrangement: {event['Location']} ({genes}{jaccard_str})\n")
                    f.write(f"        Run1: {len(event['Run1_Classifications'])} features\n")
                    f.write(f"        Run2: {len(event['Run2_Classifications'])} features\n")
                f.write("\n")

            f.write("\n")


def write_detailed_table(events: List[Dict], output_file: str):
    """
    Write detailed table in CSV format.
    """
    rows = []

    for event in events:
        # Handle multiple features by creating separate rows or combining info
        run1_feat_str = "; ".join(event['Run1_Features']) if event['Run1_Features'] else "None"
        run2_feat_str = "; ".join(event['Run2_Features']) if event['Run2_Features'] else "None"
        run1_class_str = "; ".join(event['Run1_Classifications']) if event['Run1_Classifications'] else "None"
        run2_class_str = "; ".join(event['Run2_Classifications']) if event['Run2_Classifications'] else "None"
        run1_cn_str = "; ".join([f"{cn:.1f}" for cn in event['Run1_Copy_Numbers']]) if event[
            'Run1_Copy_Numbers'] else "None"
        run2_cn_str = "; ".join([f"{cn:.1f}" for cn in event['Run2_Copy_Numbers']]) if event[
            'Run2_Copy_Numbers'] else "None"
        genes_str = "; ".join(event['Representative_Genes'])

        # Format Jaccard overlap
        jaccard_str = f"{event['Jaccard_Overlap']:.3f}" if event['Jaccard_Overlap'] is not None else "N/A"

        rows.append({
            'Sample': event['Sample'],
            'Event_Type': event['Event_Type'],
            'Run1_Features': run1_feat_str,
            'Run2_Features': run2_feat_str,
            'Run1_Classifications': run1_class_str,
            'Run2_Classifications': run2_class_str,
            'Jaccard_Overlap': jaccard_str,
            'Run1_Copy_Numbers': run1_cn_str,
            'Run2_Copy_Numbers': run2_cn_str,
            'Representative_Genes': genes_str,
            'Location': event['Location'],
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)


def generate_classification_summary(events: List[Dict]) -> Dict[str, Dict[str, int]]:
    """
    Generate summary of how each classification type changed between runs.
    Returns nested dict: {classification: {event_type: count}}
    """
    summary = defaultdict(lambda: defaultdict(int))

    for event in events:
        event_type = event['Event_Type']

        # For each run1 classification, count what happened to it
        for run1_class in event['Run1_Classifications']:
            if run1_class and run1_class != 'None':
                summary[run1_class][event_type] += 1

    return dict(summary)


def write_classification_summary(events: List[Dict], output_file: str):
    """
    Write classification change summary to file.
    """
    summary = generate_classification_summary(events)

    with open(output_file, 'w') as f:
        f.write("Classification Change Summary\n")
        f.write("=" * 40 + "\n\n")

        if not summary:
            f.write("No classification changes found.\n")
            return

        for classification in sorted(summary.keys()):
            f.write(f"Run1 {classification}:\n")

            events_for_class = summary[classification]
            total = sum(events_for_class.values())

            for event_type in ['STABLE', 'CLASS_CHANGE', 'SPLIT', 'MERGE', 'LOSS', 'COMPLEX_MERGE_SPLIT']:
                count = events_for_class.get(event_type, 0)
                if count > 0:
                    f.write(f"  {count} {event_type.lower()}\n")

            f.write(f"  Total: {total}\n\n")


def main():
    parser = argparse.ArgumentParser(description="Compare two AmpliconClassifier results tables")
    parser.add_argument("run1_file", help="Path to first AmpliconClassifier results table")
    parser.add_argument("run2_file", help="Path to second AmpliconClassifier results table")
    parser.add_argument("-o", "--output-prefix", default="comparison",
                        help="Output filename prefix (default: comparison)")

    args = parser.parse_args()

    try:
        # Read input files
        print(f"Reading {args.run1_file}...")
        df1 = pd.read_csv(args.run1_file, sep='\t')
        print(f"Reading {args.run2_file}...")
        df2 = pd.read_csv(args.run2_file, sep='\t')

        # Filter out invalid features before reporting counts
        valid_df1 = df1[df1.apply(lambda row: is_valid_feature(row, 1), axis=1)]
        valid_df2 = df2[df2.apply(lambda row: is_valid_feature(row, 2), axis=1)]

        print(f"Run 1: {len(valid_df1)} valid features across {valid_df1['Sample name'].nunique()} samples")
        print(f"Run 2: {len(valid_df2)} valid features across {valid_df2['Sample name'].nunique()} samples")
        print(f"Run 1 total entries (including no-feature samples): {len(df1)}")
        print(f"Run 2 total entries (including no-feature samples): {len(df2)}")

        # Create bipartite mapping
        print("Creating bipartite mapping...")
        mapping = create_bipartite_mapping(df1, df2)

        # Analyze events
        print("Analyzing events...")
        events = analyze_events(mapping)

        print(f"Found {len(events)} events")

        # Write outputs
        report_file = f"{args.output_prefix}_comparison_report.txt"
        table_file = f"{args.output_prefix}_comparison_table.csv"

        print(f"Writing hierarchical report to {report_file}...")
        write_hierarchical_report(events, report_file)

        print(f"Writing detailed table to {table_file}...")
        write_detailed_table(events, table_file)

        classification_summary_file = f"{args.output_prefix}_classification_summary.txt"
        print(f"Writing classification summary to {classification_summary_file}...")
        write_classification_summary(events, classification_summary_file)

        # Print summary
        event_counts = Counter(event['Event_Type'] for event in events)

        change_types = ['GAIN', 'LOSS', 'CLASS_CHANGE', 'SPLIT', 'MERGE', 'COMPLEX_MERGE_SPLIT']
        change_events = [e for e in events if e['Event_Type'] in change_types]
        stable_count = event_counts.get('STABLE', 0)

        samples_with_changes = len(set(event['Sample'] for event in change_events))
        samples_with_stable = len(set(event['Sample'] for event in events if event['Event_Type'] == 'STABLE'))

        print(f"\nSummary:")
        if change_events:
            print(f"  Changes: {len(change_events)} events across {samples_with_changes} samples")
            for event_type in change_types:
                count = event_counts.get(event_type, 0)
                if count > 0:
                    print(f"    {count} {event_type.lower()}")
        if stable_count > 0:
            print(f"  Stable: {stable_count} features across {samples_with_stable} samples")

        # Check for unmatched samples
        samples1 = set(df1['Sample name'].unique())
        samples2 = set(df2['Sample name'].unique())
        unmatched1 = samples1 - samples2
        unmatched2 = samples2 - samples1

        if unmatched1:
            print(f"Samples only in run 1: {', '.join(sorted(unmatched1))}")
        if unmatched2:
            print(f"Samples only in run 2: {', '.join(sorted(unmatched2))}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()