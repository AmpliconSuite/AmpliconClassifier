#!/usr/bin/env python3
"""
Organize amplicon images by classification into directories with symlinks or copies.

This script reads amplicon classification data and creates organized directories
containing symlinks (or copies) of PNG/PDF files based on their classifications.
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

import pandas as pd


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Organize amplicon images by classification",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.txt classifications.tsv --classes ecdna bfb
  %(prog)s input.txt classifications.tsv --classes all --copy
  %(prog)s input.txt classifications.tsv --classes ecdna linear --output my_results/
        """
    )

    parser.add_argument(
        'input_file',
        help='Input file mapping sample names to data locations'
    )

    parser.add_argument(
        'classification_file',
        help='TSV file with amplicon classifications (amplicon_classification_profiles.tsv)'
    )

    parser.add_argument(
        '--classes',
        nargs='+',
        default=['ecdna', 'bfb'],
        help='Classification categories to organize (default: ecdna bfb). '
             'Options: ecdna, bfb, complex-non-cyclic, linear, no-amp-invalid, all'
    )

    parser.add_argument(
        '--copy',
        action='store_true',
        help='Copy files instead of creating symlinks (default: create symlinks)'
    )

    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output directory (default: derived from input filename)'
    )

    parser.add_argument(
        '--file-types',
        nargs='+',
        default=['png'],
        choices=['png', 'pdf'],
        help='File types to organize (default: png)'
    )

    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite output directory if it exists'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print verbose output'
    )

    return parser.parse_args()


def parse_input_file(input_file, verbose=False):
    """
    Parse the input file to extract sample names and their data locations.

    Args:
        input_file: Path to input file
        verbose: Print verbose output

    Returns:
        Dictionary mapping sample names to their data directory paths
    """
    data_loc_dict = {}

    if not os.path.exists(input_file):
        print(f"ERROR: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(f"Reading input file: {input_file}")

    with open(input_file) as infile:
        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 2:
                print(f"WARNING: Line {line_num} has fewer than 2 fields, skipping", file=sys.stderr)
                continue

            # Extract sample name (remove _amplicon suffix if present)
            sample_name = fields[0]
            if "_amplicon" in sample_name:
                sample_name = sample_name.split("_amplicon")[0]

            # Extract directory path from second field
            file_path = fields[1]
            data_dir = str(Path(file_path).parent) + "/"

            data_loc_dict[sample_name] = data_dir

    if verbose:
        print(f"Found {len(data_loc_dict)} unique samples")

    return data_loc_dict


def normalize_class_name(class_name):
    """
    Normalize classification names to consistent format.

    Args:
        class_name: Classification name from user input or data

    Returns:
        Normalized classification name
    """
    normalize_map = {
        'ecdna': 'ecdna',
        'bfb': 'bfb',
        'complex-non-cyclic': 'complex-non-cyclic',
        'cnc': 'complex-non-cyclic',
        'linear': 'linear',
        'no-amp-invalid': 'no-amp-invalid',
        'no amp/invalid': 'no-amp-invalid',
        'invalid': 'no-amp-invalid',
    }

    return normalize_map.get(class_name.lower(), class_name.lower())


def create_output_directories(output_dir, classes, file_types, force=False, verbose=False):
    """
    Create output directory structure.

    Args:
        output_dir: Base output directory
        classes: List of classification categories
        file_types: List of file types (png, pdf)
        force: Overwrite if exists
        verbose: Print verbose output
    """
    if os.path.exists(output_dir):
        if not force:
            print(f"ERROR: Output directory '{output_dir}' already exists. Use --force to overwrite.",
                  file=sys.stderr)
            sys.exit(1)
        else:
            if verbose:
                print(f"Removing existing directory: {output_dir}")
            shutil.rmtree(output_dir)

    if verbose:
        print(f"Creating output directory structure in: {output_dir}")

    os.makedirs(output_dir)

    for class_name in classes:
        for file_type in file_types:
            subdir = os.path.join(output_dir, class_name, f"{file_type}_links")
            os.makedirs(subdir)
            if verbose:
                print(f"  Created: {subdir}")


def get_classification_category(row):
    """
    Determine classification category for an amplicon.

    Priority:
    1. ecDNA+ if Positive
    2. BFB+ if Positive
    3. Otherwise use amplicon_decomposition_class

    Args:
        row: DataFrame row with classification data

    Returns:
        List of classification categories (can be multiple)
    """
    categories = []

    # Check ecDNA+ and BFB+ status
    if row.get("ecDNA+") == "Positive":
        categories.append("ecdna")

    if row.get("BFB+") == "Positive":
        categories.append("bfb")

    # If neither ecDNA+ nor BFB+ is positive, use amplicon_decomposition_class
    if not categories:
        amp_class = row.get("amplicon_decomposition_class", "")
        if amp_class:
            normalized = normalize_class_name(amp_class)
            categories.append(normalized)

    return categories


def process_classifications(classification_file, data_loc_dict, output_dir,
                            requested_classes, file_types, use_copy=False, verbose=False):
    """
    Process classification file and create symlinks/copies.

    Args:
        classification_file: Path to TSV with classifications
        data_loc_dict: Dictionary of sample names to data locations
        output_dir: Base output directory
        requested_classes: List of classes to process
        file_types: List of file types to process
        use_copy: If True, copy files; if False, create symlinks
        verbose: Print verbose output
    """
    if not os.path.exists(classification_file):
        print(f"ERROR: Classification file '{classification_file}' not found", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(f"Reading classification file: {classification_file}")

    # Read classification data
    try:
        class_data = pd.read_csv(classification_file, sep="\t")
    except Exception as e:
        print(f"ERROR: Failed to read classification file: {e}", file=sys.stderr)
        sys.exit(1)

    # Verify required columns
    required_cols = ["sample_name", "amplicon_number"]
    missing_cols = [col for col in required_cols if col not in class_data.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)

    # Normalize requested classes
    normalized_requested = set()
    if 'all' in [c.lower() for c in requested_classes]:
        normalized_requested = {'ecdna', 'bfb', 'complex-non-cyclic', 'linear', 'no-amp-invalid'}
    else:
        normalized_requested = {normalize_class_name(c) for c in requested_classes}

    # Statistics
    stats = {
        'total': 0,
        'processed': 0,
        'skipped_no_match': 0,
        'missing_files': 0,
        'errors': 0
    }

    action_word = "Copying" if use_copy else "Linking"

    # Process each amplicon
    for idx, row in class_data.iterrows():
        stats['total'] += 1

        # Extract sample info
        sample_name = str(row["sample_name"])
        if "_amplicon" in sample_name:
            sample_name = sample_name.split("_amplicon")[0]

        amplicon_num = str(row["amplicon_number"])

        # Get classification categories
        categories = get_classification_category(row)

        # Filter to requested classes
        categories_to_process = [c for c in categories if c in normalized_requested]

        if not categories_to_process:
            stats['skipped_no_match'] += 1
            continue

        # Get data location
        if sample_name not in data_loc_dict:
            print(f"WARNING: Sample '{sample_name}' not found in input file, skipping",
                  file=sys.stderr)
            stats['errors'] += 1
            continue

        data_dir = data_loc_dict[sample_name]

        # Process each file type
        for file_type in file_types:
            filename = f"{sample_name}_{amplicon_num}.{file_type}"
            source_path = os.path.join(data_dir, filename)

            # Check if source file exists
            if not os.path.exists(source_path):
                print(f"WARNING: Source file not found: {source_path}", file=sys.stderr)
                stats['missing_files'] += 1
                continue

            # Create link/copy for each applicable category
            for category in categories_to_process:
                dest_dir = os.path.join(output_dir, category, f"{file_type}_links")
                dest_path = os.path.join(dest_dir, filename)

                try:
                    if use_copy:
                        shutil.copy2(source_path, dest_path)
                    else:
                        # Create relative or absolute symlink
                        os.symlink(os.path.abspath(source_path), dest_path)

                    if verbose:
                        print(f"  {action_word}: {filename} -> {category}/{file_type}_links/")

                    stats['processed'] += 1

                except Exception as e:
                    print(f"ERROR: Failed to process {source_path} -> {dest_path}: {e}",
                          file=sys.stderr)
                    stats['errors'] += 1

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total amplicons in classification file: {stats['total']}")
    print(f"Files successfully processed: {stats['processed']}")
    print(f"Amplicons skipped (not in requested classes): {stats['skipped_no_match']}")
    print(f"Missing source files: {stats['missing_files']}")
    print(f"Errors: {stats['errors']}")
    print(f"\nOutput directory: {output_dir}")


def main():
    """Main function."""
    args = parse_arguments()

    # Determine output directory
    if args.output:
        output_dir = args.output
    else:
        # Derive from input filename
        base_name = Path(args.input_file).stem
        output_dir = f"{base_name}_images"

    if args.verbose:
        print("=" * 60)
        print("AMPLICON IMAGE ORGANIZER")
        print("=" * 60)
        print(f"Input file: {args.input_file}")
        print(f"Classification file: {args.classification_file}")
        print(f"Output directory: {output_dir}")
        print(f"Classes to organize: {', '.join(args.classes)}")
        print(f"File types: {', '.join(args.file_types)}")
        print(f"Mode: {'Copy' if args.copy else 'Symlink'}")
        print()

    # Parse input file
    data_loc_dict = parse_input_file(args.input_file, verbose=args.verbose)

    # Normalize class names for directory creation
    normalized_classes = set()
    if 'all' in [c.lower() for c in args.classes]:
        normalized_classes = {'ecdna', 'bfb', 'complex-non-cyclic', 'linear', 'no-amp-invalid'}
    else:
        normalized_classes = {normalize_class_name(c) for c in args.classes}

    # Create output directory structure
    create_output_directories(output_dir, normalized_classes, args.file_types,
                              force=args.force, verbose=args.verbose)

    # Process classifications and create links/copies
    process_classifications(args.classification_file, data_loc_dict, output_dir,
                            args.classes, args.file_types, use_copy=args.copy,
                            verbose=args.verbose)


if __name__ == "__main__":
    main()