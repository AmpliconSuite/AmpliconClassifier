#!/usr/bin/env python3

"""
Shared utilities for extracting and analyzing read support evidence from BAM files.
Used by profile_connections.py and validate_sv_calls.py.

Author: Jens Luebeck (jluebeck [at] ucsd.edu)
"""

import re
from collections import defaultdict

import pysam
from intervaltree import IntervalTree

# Global tolerance for breakpoint matching (bp)
bp_tol = 150


# =============================================================================
# Chromosome and Coordinate Utilities
# =============================================================================

def get_chrom_value(chrom):
    """
    Convert chromosome name to a sortable numeric value.

    Args:
        chrom (str): Chromosome name (with or without 'chr' prefix)

    Returns:
        int: Numeric value for sorting chromosomes
    """
    # Strip 'chr' prefix if present
    if chrom.startswith("chr"):
        val = chrom[3:]
    else:
        val = chrom
    # Convert to comparable value
    return int(val) if val.isnumeric() else ord(val)


def is_valid_chrom(chrom):
    """
    Check if chromosome name is valid (excludes mitochondrial DNA).

    Args:
        chrom (str): Chromosome name

    Returns:
        bool: True if valid chromosome, False otherwise
    """
    # TODO: Enable functionality with viral samples.

    if chrom.startswith("chr"):
        val = chrom[3:]
    else:
        val = chrom
    # Convert to comparable value
    if val.isnumeric():
        return int(val) > 0

    else:
        if val == "M" or val == "MT":
            return False

        try:
            letter = ord(val)
            return True

        except TypeError:
            return False


def is_coordinate_less_than(chrom1, pos1, chrom2, pos2):
    """
    Compare if first genome coordinate is less than second genome coordinate.

    Args:
        chrom1 (str): First chromosome name
        pos1 (int): First position
        chrom2 (str): Second chromosome name
        pos2 (int): Second position

    Returns:
        bool: True if (chrom1, pos1) < (chrom2, pos2), False otherwise
    """
    # Compare chromosome values using shared function
    chrom1_val = get_chrom_value(chrom1)
    chrom2_val = get_chrom_value(chrom2)

    if chrom1_val < chrom2_val:
        return True
    elif chrom1_val > chrom2_val:
        return False
    else:
        # If chromosomes are the same, compare positions
        return pos1 < pos2


def sort_coordinate_tuples(tuple_list):
    """
    Sort list of tuples containing two genome coordinates and additional fields.

    Args:
        tuple_list (list): List of tuples in format:
            [(chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type), ...]

    Returns:
        list: Sorted list of tuples
    """

    def tuple_sort_key(t):
        chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, left_mapq, right_mapq = t

        # Use shared function to convert chromosome names
        chrom1_val = get_chrom_value(chrom1)
        chrom2_val = get_chrom_value(chrom2)
        if left_mapq is None:
            left_mapq = -1
        if right_mapq is None:
            right_mapq = -1

        # Return composite key with converted chromosome values
        return chrom1_val, chrom2_val, pos1, pos2, dir1, dir2, ev_type, left_mapq, right_mapq

    return sorted(tuple_list, key=tuple_sort_key)


# =============================================================================
# CIGAR and Alignment Processing
# =============================================================================

def parse_sa_cigar(cigar_str):
    """
    Parse CIGAR string from SA tag into tuples.

    Args:
        cigar_str (str): CIGAR string from SA tag

    Returns:
        list: List of (operation, length) tuples
    """
    cigar_tuples = []
    for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar_str):
        length = int(match.group(1))
        op = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5,
              'P': 6, '=': 7, 'X': 8}[match.group(2)]
        cigar_tuples.append((op, length))
    return cigar_tuples


def get_query_positions(cigar, is_primary, read, is_reverse=None, ref_start=None):
    """
    Calculate reference and query coordinates for an alignment.

    Args:
        cigar: List of CIGAR tuples (op, length)
        is_primary: Boolean indicating if this is the primary alignment
        read: pysam alignment object (used for primary alignment info)
        is_reverse: Boolean indicating if alignment is reverse complemented (for supplementary)
        ref_start: Reference start position (for supplementary)

    Returns:
        tuple: (query_start, query_end, ref_start, ref_end)
    """
    # Use pysam properties for primary alignment reference coordinates
    if is_primary:
        ref_start = read.reference_start
        ref_end = read.reference_end
        is_reverse = read.is_reverse
    else:
        ref_end = ref_start
        for op, length in cigar:
            if op in [0, 2, 3, 7, 8]:  # M, D, N, =, X all advance ref position
                ref_end += length

    # Calculate total query length including hard clips
    total_query_length = sum(length for op, length in cigar if op in [0, 1, 4, 5])  # M, I, S, H

    # Find start of alignment on query (after leading clips)
    query_start = 0
    for op, length in cigar:
        if op in [4, 5]:  # S, H
            query_start += length
        else:
            break

    # Find end of alignment on query (before trailing clips)
    query_end = total_query_length
    for op, length in reversed(cigar):
        if op in [4, 5]:  # S, H
            query_end -= length
        else:
            break

    # Reverse positions if alignment is reverse complemented
    if is_reverse:
        query_start, query_end = total_query_length - query_end, total_query_length - query_start

    return query_start, query_end, ref_start, ref_end


def calculate_query_overlap(primary_qstart, primary_qend, supp_qstart, supp_qend):
    """
    Calculate query overlap between primary and supplementary alignments.

    Args:
        primary_qstart (int): Primary alignment query start
        primary_qend (int): Primary alignment query end
        supp_qstart (int): Supplementary alignment query start
        supp_qend (int): Supplementary alignment query end

    Returns:
        int: Number of overlapping bases in query coordinates
    """
    # Query coordinates are now handled in get_query_positions
    # Just calculate direct overlap
    overlap_start = max(primary_qstart, supp_qstart)
    overlap_end = min(primary_qend, supp_qend)
    return max(0, overlap_end - overlap_start)


# =============================================================================
# Split Read Evidence Extraction
# =============================================================================

def find_split_read_evidence(read, lcD, min_sv_size=5000, min_aln_length=20, max_overlap=10, min_mapq=5, debug=False):
    """
    Identify split read evidence for SVs from a single read.

    Args:
        read: pysam alignment object of primary alignment
        lcD: low complexity region dictionary
        min_sv_size: minimum size for intrachromosomal SVs
        min_aln_length: minimum length of clip to consider
        max_overlap: maximum allowed overlap between alignments
        min_mapq: minimum mapping quality
        debug: whether to print debug information

    Returns:
        tuple or None: (chrom1, pos1, chrom2, pos2, dir1, dir2, 'split', mapq1, mapq2) or None
    """
    evidence = None
    if not read.cigartuples:
        return evidence

    # Only process primary alignments with SA tags
    if not read.has_tag('SA'):
        return evidence

    # Calculate total query length including hard clips
    total_query_length = sum(l for op, l in read.cigartuples if op in [0, 1, 4, 5])  # M, I, S, H

    sa_info = read.get_tag('SA').split(',')
    primary_chrom = read.reference_name
    supp_chrom = sa_info[0]
    supp_pos = int(sa_info[1])
    supp_is_reverse = sa_info[2] == '-'
    supp_cigar = parse_sa_cigar(sa_info[3])
    supp_mapq = int(sa_info[4])

    if not is_valid_chrom(primary_chrom) or not is_valid_chrom(supp_chrom) or not supp_mapq > min_mapq:
        return evidence

    # Store coordinates for both alignments
    primary_qstart, primary_qend, primary_rstart, primary_rend = get_query_positions(
        read.cigartuples, True, read)
    supp_qstart, supp_qend, supp_rstart, supp_rend = get_query_positions(
        supp_cigar, False, read, supp_is_reverse, supp_pos)

    # Too short
    if supp_qend - supp_qstart < min_aln_length:
        return evidence

    if debug:
        print("\nProcessing read:", read.query_name)
        print("Is read1:", read.is_read1)
        print("Primary alignment:", read.reference_name, read.reference_start, read.reference_end,
              "reverse" if read.is_reverse else "forward")
        print("Primary CIGAR:", read.cigarstring)
        print("Supplementary CIGAR:", supp_cigar)
        print("Total query length (including hard clips):", total_query_length)

        print("\nAlignment coordinates:")
        print("Primary alignment:")
        print("  Reference: {}:{}-{}".format(primary_chrom, primary_rstart, primary_rend))
        print("  Query mapping:")
        print("    ref:{} -> query:{}".format(primary_rstart, primary_qstart))
        print("    ref:{} -> query:{}".format(primary_rend, primary_qend))
        print("  Orientation:", "reverse" if read.is_reverse else "forward")

        print("\nSupplementary alignment:")
        print("  Reference: {}:{}-{}".format(supp_chrom, supp_rstart, supp_rend))
        print("  Query mapping:")
        print("    ref:{} -> query:{}".format(supp_rstart, supp_qstart))
        print("    ref:{} -> query:{}".format(supp_rend, supp_qend))
        print("  Orientation:", "reverse" if supp_is_reverse else "forward")

    # Calculate overlap
    query_overlap = calculate_query_overlap(primary_qstart, primary_qend, supp_qstart, supp_qend)

    if query_overlap > max_overlap:
        if debug:
            print("Rejecting: Query overlap too large:", query_overlap)
        return evidence

    dir1 = '-' if read.is_reverse else '+'
    dir2 = '-' if supp_is_reverse else '+'

    primary_mapq, supp_mapq = read.mapq, supp_mapq

    # if supp comes before primary (left clip) swap primary and supplementary so they are ordered on the molecule
    if supp_qstart < primary_qstart:
        temp_qstart, temp_qend, temp_rstart, temp_rend = primary_qstart, primary_qend, primary_rstart, primary_rend
        primary_qstart, primary_qend, primary_rstart, primary_rend = supp_qstart, supp_qend, supp_rstart, supp_rend
        supp_qstart, supp_qend, supp_rstart, supp_rend = temp_qstart, temp_qend, temp_rstart, temp_rend
        dir1, dir2 = dir2, dir1
        primary_chrom, supp_chrom = supp_chrom, primary_chrom
        primary_mapq, supp_mapq = supp_mapq, primary_mapq

    if dir1 == "+":
        primary_break = primary_rend
    else:
        primary_break = primary_rstart

    if dir2 == "+":
        supp_break = supp_rstart
    else:
        supp_break = supp_rend

    if not (lcD[primary_chrom][primary_rstart:primary_rend] or lcD[supp_chrom][supp_rstart:supp_rend]):
        if primary_chrom == supp_chrom and abs(primary_break - supp_break) < min_sv_size and dir1 == dir2:
            if debug:
                print("  Rejected: Same chromosome and too close")
            return evidence

        if debug:
            print("  Added evidence: ")
            print("   " + str((primary_chrom, primary_break, supp_chrom, supp_break, dir1, dir2, 'split', primary_mapq,
                               supp_mapq)))
        evidence = (primary_chrom, primary_break, supp_chrom, supp_break, dir1, dir2, 'split', primary_mapq, supp_mapq)

    return evidence


# =============================================================================
# Discordant Pair Evidence Extraction
# =============================================================================

def create_merged_mapq_windows(mate_locations, window_size=50, extend_size=10000, verbose=False):
    """
    Create merged genomic windows from mate locations.
    Single reads get small windows of size window_size,
    while clusters of reads within extend_size get merged into larger windows.
    The final position in each cluster gets minimal extension.

    Args:
        mate_locations: List of (chrom, pos, query_name, is_read2) tuples
        window_size: Size around isolated positions
        extend_size: Size to extend around each position for merging

    Returns:
        List of (chrom, start, end) tuples representing merged windows
    """
    # Group by chromosome
    chrom_groups = {}
    for chrom, pos, _, _ in mate_locations:
        if chrom not in chrom_groups:
            chrom_groups[chrom] = []
        chrom_groups[chrom].append(pos)

    merged_windows = []
    # Process each chromosome separately
    for chrom, positions in chrom_groups.items():
        if not positions:
            continue

        # Sort positions
        positions.sort()

        # Initialize first window
        curr_start = positions[0]
        curr_end = positions[0] + window_size
        last_pos = positions[0]

        # Merge windows
        for pos in positions[1:]:
            if pos <= curr_end + extend_size:
                # Update last position and extend window
                last_pos = pos
                curr_end = max(curr_end, pos + window_size)
            else:
                # Save current window with minimal extension after last position
                merged_windows.append((chrom, max(0, curr_start - 1), last_pos + 1))
                curr_start = pos
                curr_end = pos + window_size
                last_pos = pos

        # Add final window with minimal extension
        merged_windows.append((chrom, max(0, curr_start - 1), last_pos + 1))

    if verbose:
        print("Windows to check for mates: {}".format(len(merged_windows)))

    return merged_windows


def get_mate_mapqs(bam_path, mate_locations):
    """
    Efficiently fetch mapping qualities for all mate reads.

    Args:
        bam_path: Path to BAM file
        mate_locations: List of (chrom, pos, query_name, is_read2) tuples

    Returns:
        dict: Mapping of query_names to mate mapqs
    """
    # Create merged windows for efficient fetching
    windows = create_merged_mapq_windows(mate_locations)

    # Create mapping of read names to their mate info for matching
    mate_map = {name: (chrom, pos, is_read2)
                for chrom, pos, name, is_read2 in mate_locations}

    # Dictionary to store mate mapqs
    mate_mapqs = {}

    bamfile = pysam.AlignmentFile(bam_path, "rb")

    # Process each merged window
    for chrom, start, end in windows:
        for read in bamfile.fetch(chrom, start, end):
            # Check if this read is one of our mates
            if read.query_name in mate_map:
                mate_info = mate_map[read.query_name]
                # Verify this is the correct mate read
                if (read.reference_name == mate_info[0] and
                        read.reference_start == mate_info[1] and
                        read.is_read2 == mate_info[2]):
                    mate_mapqs[read.query_name] = read.mapping_quality

    bamfile.close()
    return mate_mapqs


def find_discordant_pair_evidence(read, min_sv_size, lcD, debug=False):
    """
    Identify discordant pair evidence for SVs from a single read.
    Uses proper breakpoint positions based on orientation.

    Args:
        read: pysam alignment object
        min_sv_size: minimum size for intrachromosomal SVs
        lcD: low complexity region dictionary
        debug: whether to print debug information

    Returns:
        tuple: (evidence, mate_info) where evidence is a list or None,
               and mate_info is (chrom, pos, query_name, is_read2) or None
    """
    if not (read.is_paired and read.mate_is_mapped and
            not read.is_secondary and not read.is_supplementary and
            read.query_name):
        return None, None

    # Get basic information about the read pair
    read1_chrom = read.reference_name
    read2_chrom = read.next_reference_name

    if not is_valid_chrom(read1_chrom) or not is_valid_chrom(read2_chrom):
        return None, None

    # Get proper positions based on orientation
    if not read.is_reverse:
        init_pos1 = read.reference_end  # If read1 is forward, use end
    else:
        init_pos1 = read.reference_start  # If read1 is reverse, use end

    if read.mate_is_reverse:
        init_pos2 = read.next_reference_start  # If read2 is reverse, use start
    else:
        # For mate end position, we need to estimate since we don't have the actual alignment
        init_pos2 = read.next_reference_start + read.query_length

    if lcD[read1_chrom][read.reference_start:read.reference_end] or \
            lcD[read2_chrom][read.next_reference_start:read.next_reference_start + read.query_length]:
        return None, None

    # Determine orientations
    dir1 = '-' if read.is_reverse else '+'
    dir2 = '-' if read.mate_is_reverse else '+'

    mapq1, mapq2 = read.mapq, None

    # re-order read1, read2 so that read1 is always on the smaller coordinate.
    swap = False
    if is_coordinate_less_than(read2_chrom, init_pos2, read1_chrom, init_pos1):
        swap = True
        read1_chrom, read2_chrom = read2_chrom, read1_chrom
        dir1, dir2 = dir2, dir1
        mapq1, mapq2 = mapq2, mapq1

    # convert orientation to a genome flow direction
    if dir1 == "+" and dir2 == "-":
        dir1, dir2 = "+", "+"

    elif dir1 == "+" and dir2 == "+":
        dir1, dir2 = "+", "-"

    elif dir1 == "-" and dir2 == "-":
        dir1, dir2 = "-", "+"

    elif dir1 == "-" and dir2 == "+":
        dir1, dir2 = "-", "-"

    # now figure out which start/end you need
    if not swap:
        if dir1 == "+":
            pos1 = read.reference_end
        else:
            pos1 = read.reference_start

        if dir2 == "+":  # flow orientation, not illumina orientation
            pos2 = read.next_reference_start
        else:
            pos2 = read.next_reference_start + read.query_length

    else:
        if dir1 == "+":
            pos1 = read.next_reference_start + read.query_length

        else:
            pos1 = read.next_reference_start

        if dir2 == "+":
            pos2 = read.reference_start
        else:
            pos2 = read.reference_end

    # Different chromosomes or same chromosome with unexpected orientation/size
    if (read1_chrom != read2_chrom or
            abs(read.template_length) >= min_sv_size or
            read.is_reverse == read.mate_is_reverse):  # Check for non-FR pairs

        # populate the missing mapq
        mate_to_get = (read.next_reference_name, read.next_reference_start, read.query_name, not read.is_read2)
        return [read1_chrom, pos1, read2_chrom, pos2, dir1, dir2, 'discordant', mapq1, mapq2], mate_to_get

    return None, None


# =============================================================================
# SV Matching Utilities
# =============================================================================

def build_sv_interval_dict(sv_list, window=bp_tol):
    """
    Build interval tree dictionary for efficient SV lookup.

    Args:
        sv_list: List of SV objects with lchrom, lpos, rchrom, rpos attributes
        window: Window size for interval matching (default: bp_tol)

    Returns:
        defaultdict: Maps chromosome to IntervalTree containing SVs
    """
    sv_ivald = defaultdict(IntervalTree)
    for sv in sv_list:
        sv_ivald[sv.lchrom].addi(sv.lpos - window, sv.lpos + window, sv)
        sv_ivald[sv.rchrom].addi(sv.rpos - window, sv.rpos + window, sv)

    return sv_ivald


def match_evidence_to_sv(evidence, sv_ivald, window=bp_tol):
    """
    Match read evidence to an SV in the interval dictionary.

    Args:
        evidence: Tuple of (chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, mapq1, mapq2)
        sv_ivald: Interval dictionary from build_sv_interval_dict
        window: Tolerance window for matching (default: bp_tol)

    Returns:
        Matched SV object or None if no match found
    """
    chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, mapq1, mapq2 = evidence

    # Find potential matches on both sides
    l_hits = sv_ivald[chrom1][pos1]
    r_hits = sv_ivald[chrom2][pos2]
    l_r_hits = set(x.data for x in l_hits) & set(x.data for x in r_hits)

    for sv in l_r_hits:
        # Make sure it's not hitting the same side twice
        lbp, rbp = sorted(((chrom1, pos1), (chrom2, pos2)))
        sv_lbp, sv_rbp = sorted(((sv.lchrom, sv.lpos), (sv.rchrom, sv.rpos)))

        if lbp[0] == sv_lbp[0] and rbp[0] == sv_rbp[0]:
            if abs(lbp[1] - sv_lbp[1]) <= window and abs(rbp[1] - sv_rbp[1]) <= window:
                # TODO: Add directionality checks here later
                return sv

    return None


# =============================================================================
# Main Read Extraction Function
# =============================================================================

def find_sv_supporting_reads(bam_path, chrom, start, end, lcD, min_mapq=5, min_sv_size=5000, verbose=False,
                             include_read_names=False):
    """
    Find all split read and discordant pair evidence in a genomic region.

    Args:
        bam_path: Path to BAM file
        chrom: Chromosome name
        start: Start position
        end: End position
        lcD: Low complexity region dictionary
        min_mapq: Minimum mapping quality
        min_sv_size: Minimum SV size for discordant pairs
        verbose: Print verbose output
        include_read_names: If True, append read_name to evidence tuples

    Returns:
        list: Combined list of split read and discordant pair evidence tuples
              Format: (chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, mapq1, mapq2[, read_name])
              read_name is appended only if include_read_names=True
    """
    split_evidence = []
    disc_evidence = []
    mate_locations = []  # List of (chrom, pos, query_name, is_read2) tuples
    used_reads = set()

    bamfile = pysam.AlignmentFile(bam_path, "rb")
    for read in bamfile.fetch(chrom, start - bp_tol, end + bp_tol):
        if read.mapping_quality < min_mapq:
            continue

        # Split reads
        split_ev = find_split_read_evidence(read, lcD, debug=verbose)
        if split_ev:
            if include_read_names:
                split_ev = split_ev + (read.query_name,)
            split_evidence.append(split_ev)

        # Discordant pairs
        if read.query_name not in used_reads:
            disc_ev, mate_info = find_discordant_pair_evidence(read, min_sv_size, lcD, debug=verbose)
            if disc_ev:
                used_reads.add(read.query_name)
                if include_read_names:
                    disc_ev.append(read.query_name)
                disc_evidence.append(disc_ev)
                mate_locations.append(mate_info)

    bamfile.close()

    fixed_mapqs = 0
    unfixed_mapqs = 0
    mate_mapq_dict = get_mate_mapqs(bam_path, mate_locations)

    updated_disc_evidence = []
    for disc_ev, mate_info in zip(disc_evidence, mate_locations):
        read_name = mate_info[2]
        mate_mapq = mate_mapq_dict.get(read_name, None)  # Default to None if not found
        if mate_mapq is not None:
            fixed_mapqs += 1
        else:
            unfixed_mapqs += 1

        # Convert to list for modification
        disc_ev_list = list(disc_ev)

        # Update MAPQ values (they're at indices 7 and 8, regardless of whether read_name is present)
        if disc_ev_list[7] is None:
            disc_ev_list[7] = mate_mapq
        else:
            disc_ev_list[8] = mate_mapq

        updated_disc_evidence.append(tuple(disc_ev_list))

    # Return combined evidence
    if verbose:
        print("Retrieved {} mate mapqs for {} discordant pairs".format(fixed_mapqs, fixed_mapqs + unfixed_mapqs))

    return split_evidence + updated_disc_evidence