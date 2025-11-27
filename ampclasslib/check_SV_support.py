#!/usr/bin/env python3

"""
Validate structural variant calls from a table by checking for supporting reads in a BAM file.

Takes an SV calls table (TSV format) and a BAM file, and annotates each SV with the number
of supporting discordant pairs and split reads found in the BAM.

Author: Jens Luebeck (jluebeck [at] ucsd.edu)
"""

import argparse
from collections import defaultdict
import os
import sys

import pysam

# Allow use as script or as library
if __name__ == "__main__":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from ac_util import set_lcd, breakpoint
    from read_support_utils import *

else:
    from ampclasslib.ac_util import set_lcd, breakpoint
    from ampclasslib.read_support_utils import *


class SVRecord:
    """
    Represents a structural variant record from the input table.
    Stores all original fields plus counters for validation support.
    Contains a breakpoint instance for matching.
    """

    def __init__(self, fields, header):
        """
        Initialize SVRecord from a line of the input TSV.

        Args:
            fields: List of field values from TSV line
            header: List of column names from TSV header
        """
        # Store original row for output
        self.original_row = fields

        # Parse key fields for matching
        self.chrom1 = fields[0]
        self.pos1 = int(fields[1])
        self.chrom2 = fields[2]
        self.pos2 = int(fields[3])

        # Store other fields if they exist
        if len(fields) > 4:
            self.sv_type = fields[4] if len(fields) > 4 else None
            self.read_support = fields[5] if len(fields) > 5 else None
            self.features = fields[6] if len(fields) > 6 else None
            self.orientation = fields[7] if len(fields) > 7 else None
            self.pos1_flanking = fields[8] if len(fields) > 8 else None
            self.pos2_flanking = fields[9] if len(fields) > 9 else None
            self.homology_length = fields[10] if len(fields) > 10 else None
            self.homology_sequence = fields[11] if len(fields) > 11 else None

        # Create breakpoint instance for matching
        # Parse orientation into ldir/rdir if present
        ldir, rdir = None, None
        if hasattr(self, 'orientation') and self.orientation and self.orientation != 'None':
            if len(self.orientation) == 2:
                ldir = self.orientation[0]
                rdir = self.orientation[1]

        self.bp = breakpoint(
            lchrom=self.chrom1,
            lpos=self.pos1,
            rchrom=self.chrom2,
            rpos=self.pos2,
            ldir=ldir,
            rdir=rdir,
            homlen=self.homology_length if hasattr(self, 'homology_length') else "NA",
            homseq=self.homology_sequence if hasattr(self, 'homology_sequence') else "NA"
        )

        # Initialize validation counters
        self.bam_discordant_pairs = 0
        self.bam_split_reads = 0
        self.bam_total_support = 0

    def __str__(self):
        """String representation for debugging."""
        return f"SV({self.chrom1}:{self.pos1}-{self.chrom2}:{self.pos2})"

    def calculate_total_support(self):
        """Calculate total support from discordant and split reads."""
        self.bam_total_support = self.bam_discordant_pairs + self.bam_split_reads


def parse_sv_table(input_file):
    """
    Parse the input SV calls table.

    Args:
        input_file: Path to input TSV file

    Returns:
        tuple: (header_list, list of SVRecord objects)
    """
    sv_records = []
    header = None

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')

            # First non-empty line is header
            if header is None:
                header = fields
                continue

            # Create SVRecord from data line
            sv_record = SVRecord(fields, header)
            sv_records.append(sv_record)

    if header is None:
        raise ValueError("Input file is empty or has no header")

    print(f"Parsed {len(sv_records)} SV records from {input_file}")
    return header, sv_records


def validate_svs_from_bam(bam_path, sv_records, lcD, tolerance=bp_tol,
                          min_mapq=5, min_sv_size=5000, window=500, verbose=False):
    """
    Validate SVs by finding supporting reads in BAM file using region-based approach.

    Args:
        bam_path: Path to BAM file
        sv_records: List of SVRecord objects to validate
        lcD: Low complexity region dictionary
        tolerance: Tolerance window for matching (default: bp_tol)
        min_mapq: Minimum mapping quality
        min_sv_size: Minimum SV size for discordant pairs
        window: Window size around each breakpoint to search (default: 500bp)
        verbose: Print verbose output
    """
    print(f"Validating {len(sv_records)} SVs against {bam_path}")
    print(f"Using {window}bp window around each breakpoint")

    total_split_matches = 0
    total_disc_matches = 0

    # Process each SV individually
    for idx, sv in enumerate(sv_records):
        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1}/{len(sv_records)} SVs...")

        if verbose:
            print(f"\nProcessing SV: {sv}")

        # Fetch reads from regions around both breakpoints
        try:
            evidence_1 = find_sv_supporting_reads(
                bam_path, sv.chrom1, sv.pos1 - window, sv.pos1 + window,
                lcD, min_mapq=min_mapq, min_sv_size=min_sv_size,
                verbose=verbose, include_read_names=True
            )
        except Exception as e:
            if verbose:
                print(f"Warning: Could not fetch reads for {sv.chrom1}:{sv.pos1}: {e}")
            evidence_1 = []

        try:
            evidence_2 = find_sv_supporting_reads(
                bam_path, sv.chrom2, sv.pos2 - window, sv.pos2 + window,
                lcD, min_mapq=min_mapq, min_sv_size=min_sv_size,
                verbose=verbose, include_read_names=True
            )
        except Exception as e:
            if verbose:
                print(f"Warning: Could not fetch reads for {sv.chrom2}:{sv.pos2}: {e}")
            evidence_2 = []

        # Combine evidence from both regions
        all_evidence = evidence_1 + evidence_2

        # Deduplicate by read name and count evidence matching this SV
        seen_read_names = set()

        for ev in all_evidence:
            # Evidence format: (chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, mapq1, mapq2, read_name)
            *coords_and_type, read_name = ev

            # Skip if we've already seen this read
            if read_name in seen_read_names:
                continue
            seen_read_names.add(read_name)

            # Check if this evidence matches THIS specific SV
            # coords_and_type is: [chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, mapq1, mapq2]
            position_match, orientation_match = match_evidence_to_sv_record(coords_and_type, sv, tolerance,
                                                                            verbose=verbose)

            # Only count if position matches AND (orientation matches OR no orientation available)
            if position_match and (orientation_match is True or orientation_match is None):
                ev_type = coords_and_type[6]
                if ev_type == 'split':
                    sv.bam_split_reads += 1
                    total_split_matches += 1
                    if verbose:
                        print(f"  Split read COUNTED: {read_name}")
                elif ev_type == 'discordant':
                    sv.bam_discordant_pairs += 1
                    total_disc_matches += 1
                    if verbose:
                        print(f"  Discordant pair COUNTED: {read_name}")
            elif position_match and orientation_match is False:
                if verbose:
                    print(f"  Evidence NOT COUNTED - orientation mismatch")

        # Calculate total support for this SV
        sv.calculate_total_support()

    print(f"\nValidation complete!")
    print(f"Total split read matches: {total_split_matches}")
    print(f"Total discordant pair matches: {total_disc_matches}")

    return sv_records


def match_evidence_to_sv_order(ev_chrom1, ev_pos1, ev_chrom2, ev_pos2, ev_dir1, ev_dir2,
                               sv_chrom1, sv_pos1, sv_chrom2, sv_pos2):
    """
    Reorder evidence coordinates to match the AA SV's coordinate order.
    The AA SV defines the coordinate relationship (pos1 vs pos2).
    We adjust the evidence to have the same relationship.

    Args:
        ev_chrom1, ev_pos1, ev_dir1: Evidence first breakpoint
        ev_chrom2, ev_pos2, ev_dir2: Evidence second breakpoint
        sv_chrom1, sv_pos1: SV first breakpoint
        sv_chrom2, sv_pos2: SV second breakpoint

    Returns:
        tuple: (chrom1, pos1, chrom2, pos2, dir1, dir2) reordered to match SV
    """
    # Sort both evidence coordinates to find left and right
    ev_sorted = sorted([(ev_chrom1, ev_pos1, ev_dir1), (ev_chrom2, ev_pos2, ev_dir2)],
                       key=lambda x: (x[0], x[1]))
    ev_left = ev_sorted[0]  # (chrom, pos, dir) for smaller coordinate
    ev_right = ev_sorted[1]  # (chrom, pos, dir) for larger coordinate

    # Sort SV coordinates to find left and right
    sv_sorted = sorted([(sv_chrom1, sv_pos1), (sv_chrom2, sv_pos2)],
                       key=lambda x: (x[0], x[1]))
    sv_left = sv_sorted[0]  # (chrom, pos) for smaller coordinate
    sv_right = sv_sorted[1]  # (chrom, pos) for larger coordinate

    # Now determine which evidence coordinate matches which SV coordinate
    # Match based on which SV position is closest to each evidence position
    dist_ev1_to_sv1 = abs(ev_pos1 - sv_pos1) if ev_chrom1 == sv_chrom1 else float('inf')
    dist_ev1_to_sv2 = abs(ev_pos1 - sv_pos2) if ev_chrom1 == sv_chrom2 else float('inf')
    dist_ev2_to_sv1 = abs(ev_pos2 - sv_pos1) if ev_chrom2 == sv_chrom1 else float('inf')
    dist_ev2_to_sv2 = abs(ev_pos2 - sv_pos2) if ev_chrom2 == sv_chrom2 else float('inf')

    # Check if ev1 matches sv1 and ev2 matches sv2, or vice versa
    if dist_ev1_to_sv1 + dist_ev2_to_sv2 < dist_ev1_to_sv2 + dist_ev2_to_sv1:
        # ev1 matches sv1, ev2 matches sv2 - keep order
        return ev_chrom1, ev_pos1, ev_chrom2, ev_pos2, ev_dir1, ev_dir2
    else:
        # ev1 matches sv2, ev2 matches sv1 - swap and flip directions
        new_dir1 = '-' if ev_dir2 == '+' else '+'
        new_dir2 = '-' if ev_dir1 == '+' else '+'
        return ev_chrom2, ev_pos2, ev_chrom1, ev_pos1, new_dir1, new_dir2


def illumina_to_sv_flow(dir1, dir2):
    """
    Convert AA's Illumina orientation to SV flow orientation.
    Flip the sign on the second direction.

    Args:
        dir1, dir2: Illumina orientation

    Returns:
        tuple: (dir1, dir2) in SV flow orientation
    """
    new_dir2 = '-' if dir2 == '+' else '+'
    return dir1, new_dir2


def match_evidence_to_sv_record(evidence, sv, tolerance, verbose=False):
    """
    Check if evidence matches a specific SV record.
    Reorders evidence to match AA's coordinate relationship, then compares orientations.

    Args:
        evidence: List [chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, mapq1, mapq2]
        sv: SVRecord object
        tolerance: Tolerance window for matching
        verbose: If True, print detailed orientation comparison

    Returns:
        tuple: (position_match, orientation_match)
               position_match: bool - True if positions match within tolerance
               orientation_match: bool or None - True if orientations match, None if no SV orientation
    """
    ev_chrom1, ev_pos1, ev_chrom2, ev_pos2 = evidence[0], evidence[1], evidence[2], evidence[3]
    ev_dir1, ev_dir2, ev_type = evidence[4], evidence[5], evidence[6]

    # Create a breakpoint from the evidence for position matching
    evidence_bp = breakpoint(
        lchrom=ev_chrom1,
        lpos=ev_pos1,
        rchrom=ev_chrom2,
        rpos=ev_pos2
    )

    # Use the breakpoint.d_similar() method for position matching
    is_position_match = sv.bp.d_similar(evidence_bp, 2 * tolerance)

    if not is_position_match:
        return False, None

    # Position matches! Now check orientation
    # Reorder evidence to match AA's coordinate order
    ev_c1, ev_p1, ev_c2, ev_p2, ev_d1, ev_d2 = match_evidence_to_sv_order(
        ev_chrom1, ev_pos1, ev_chrom2, ev_pos2, ev_dir1, ev_dir2,
        sv.chrom1, sv.pos1, sv.chrom2, sv.pos2
    )

    # Convert SV orientation from Illumina to SV flow
    sv_orientation = sv.orientation if hasattr(sv, 'orientation') and sv.orientation else None

    if sv_orientation and sv_orientation != 'None' and len(sv_orientation) == 2:
        sv_dir1, sv_dir2 = sv_orientation[0], sv_orientation[1]
        # Convert Illumina to SV flow
        sv_d1, sv_d2 = illumina_to_sv_flow(sv_dir1, sv_dir2)

        orientation_match = (sv_d1 == ev_d1 and sv_d2 == ev_d2)
    else:
        # No orientation available, can't check
        orientation_match = None

    if verbose:
        # Show detailed comparison
        print(f"\n  === Match Found ===")
        print(f"  SV Breakpoint:       {sv.chrom1}:{sv.pos1} | {sv.chrom2}:{sv.pos2}")
        print(f"  Evidence Breakpoint (original): {ev_chrom1}:{ev_pos1} | {ev_chrom2}:{ev_pos2}")
        print(f"  Evidence Breakpoint (reordered): {ev_c1}:{ev_p1} | {ev_c2}:{ev_p2}")
        print(f"  Position differences: Δ1={abs(sv.pos1 - ev_p1)}bp, Δ2={abs(sv.pos2 - ev_p2)}bp")
        print(f"  Evidence type: {ev_type}")

        # Show orientation comparison
        print(f"  SV orientation (Illumina):        {sv_orientation if sv_orientation else 'unknown'}")
        print(f"  Evidence orientation (original):  {ev_dir1}{ev_dir2}")
        print(f"  Evidence orientation (reordered): {ev_d1}{ev_d2}")

        if orientation_match is not None:
            print(f"  After conversion:")
            print(f"    SV (converted to SV flow):  {sv.chrom1}:{sv.pos1}({sv_d1}) -> {sv.chrom2}:{sv.pos2}({sv_d2})")
            print(f"    Evidence (reordered):       {ev_c1}:{ev_p1}({ev_d1}) -> {ev_c2}:{ev_p2}({ev_d2})")

            if orientation_match:
                print(f"  ✓ Orientations MATCH - counting this evidence")
            else:
                print(f"  ✗ Orientations DO NOT MATCH - NOT counting this evidence")
                print(f"     This may indicate a false positive or different SV interpretation")
        else:
            print(f"  ⚠ SV orientation not available - counting based on position only")

    return is_position_match, orientation_match


def write_output(output_file, header, sv_records):
    """
    Write annotated SV records to output file.

    Args:
        output_file: Path to output TSV file
        header: Original header list
        sv_records: List of SVRecord objects with validation counts
    """
    # Add new columns to header
    new_header = header + ['bam_discordant_pairs', 'bam_split_reads', 'bam_total_support']

    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(new_header) + '\n')

        # Write each SV record
        for sv in sv_records:
            # Original fields
            original_line = '\t'.join(sv.original_row)
            # New fields
            new_fields = [
                str(sv.bam_discordant_pairs),
                str(sv.bam_split_reads),
                str(sv.bam_total_support)
            ]
            f.write(original_line + '\t' + '\t'.join(new_fields) + '\n')

    print(f"Wrote annotated results to {output_file}")


def get_output_filename(input_file, bam_path, output_file=None):
    """
    Determine output filename based on input file and BAM name.

    Args:
        input_file: Input SV table path
        bam_path: BAM file path
        output_file: Optional user-specified output path

    Returns:
        str: Output file path
    """
    if output_file:
        return output_file

    # Extract base names
    input_base = os.path.basename(input_file)
    input_name = os.path.splitext(input_base)[0]

    bam_base = os.path.basename(bam_path)
    bam_name = bam_base.replace('.bam', '').replace('.BAM', '')

    # Construct output filename in current working directory
    output_name = f"{input_name}_annotated_{bam_name}.tsv"

    return output_name


def main():
    parser = argparse.ArgumentParser(
        description="Validate SV calls from a table by checking for supporting reads in a BAM file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    validate_sv_calls.py --input_table sv_calls.tsv --bam sample.bam --ref hg38

Output:
    Creates sv_calls_annotated_sample.tsv with three additional columns:
    - bam_discordant_pairs: Number of discordant read pairs supporting the SV
    - bam_split_reads: Number of split reads supporting the SV
    - bam_total_support: Total support (sum of above two)
        """
    )

    parser.add_argument("--input_table", "-i", required=True,
                        help="Input TSV file containing SV calls")
    parser.add_argument("--bam", "-b", required=True,
                        help="BAM file to validate SVs against")
    parser.add_argument("--ref", "-r", required=True,
                        choices=["hg19", "GRCh37", "GRCh38", "hg38", "mm10"],
                        help="Reference genome name")
    parser.add_argument("--output", "-o",
                        help="Output TSV file (default: auto-generated from input and BAM names)")
    parser.add_argument("--tolerance", "-t", type=int, default=bp_tol,
                        help=f"Breakpoint tolerance window in bp (default: {bp_tol})")
    parser.add_argument("--window", "-w", type=int, default=500,
                        help="Window size around each breakpoint to search for reads (default: 500bp)")
    parser.add_argument("--min_mapq", type=int, default=5,
                        help="Minimum mapping quality (default: 5)")
    parser.add_argument("--min_sv_size", type=int, default=0,
                        help="Minimum SV size for intrachromosomal events (default: 0)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Verbose output for debugging")

    args = parser.parse_args()

    # Normalize reference name
    if args.ref == "hg38":
        args.ref = "GRCh38"

    # Check for AA_DATA_REPO
    try:
        AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + args.ref + "/"
    except KeyError:
        sys.stderr.write("$AA_DATA_REPO not set. Please see AA installation instructions.\n")
        sys.exit(1)

    # Set up low complexity database
    lcD, cg5D = set_lcd(AA_DATA_REPO)
    print("Loaded low complexity database")

    # Determine output filename
    output_file = get_output_filename(args.input_table, args.bam, args.output)
    print(f"Output will be written to: {output_file}")

    # Parse input SV table
    header, sv_records = parse_sv_table(args.input_table)

    # Validate SVs from BAM using region-based approach
    sv_records = validate_svs_from_bam(
        args.bam, sv_records, lcD,
        tolerance=args.tolerance,
        min_mapq=args.min_mapq,
        min_sv_size=args.min_sv_size,
        window=args.window,
        verbose=args.verbose
    )

    # Write output
    write_output(output_file, header, sv_records)

    # Print summary statistics
    total_svs = len(sv_records)
    svs_with_support = sum(1 for sv in sv_records if sv.bam_total_support > 0)
    total_disc = sum(sv.bam_discordant_pairs for sv in sv_records)
    total_split = sum(sv.bam_split_reads for sv in sv_records)

    print("\n=== Summary ===")
    print(f"Total SVs in input: {total_svs}")
    print(f"SVs with support in BAM: {svs_with_support} ({100 * svs_with_support / total_svs:.1f}%)")
    print(f"Total discordant pairs: {total_disc}")
    print(f"Total split reads: {total_split}")
    print(f"Total support: {total_disc + total_split}")


if __name__ == "__main__":
    main()