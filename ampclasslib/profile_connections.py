#!/usr/bin/env python3

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"

import argparse
from collections import Counter, defaultdict
import os
import random
import re
import sys

import numpy as np
import pysam
from intervaltree import IntervalTree

# Allow use as script or as library
if __name__ == "__main__":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from ac_io import parse_bpg, bed_to_interval_dict
    from ac_util import set_lcd

else:
    from ampclasslib.ac_io import parse_bpg, bed_to_interval_dict
    from ampclasslib.ac_util import set_lcd


bp_tol = 150
random.seed(0)


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


def build_aa_bp_ivald(aa_bp_list, window=bp_tol):
    aa_bp_ivald = defaultdict(IntervalTree)
    for bp in aa_bp_list:
        aa_bp_ivald[bp.lchrom].addi(bp.lpos-window, bp.lpos+window, bp)
        aa_bp_ivald[bp.rchrom].addi(bp.rpos - window, bp.rpos + window, bp)

    return aa_bp_ivald


def translate_aa_bp_dir(bp):
    # bp is a breakpoint instance
    c1, c2 = bp.lchrom, bp.rchrom
    p1, p2 = bp.lpos, bp.rpos
    ldir, rdir = bp.ldir, bp.rdir
    if c1 == c2 and p1 > p2:
        p1, p2 = p2, p1
        ldir, rdir = rdir, ldir

    # TODO: FINISH IMPLEMENTING
    #new_rdir = "+" if rdir == "-" else "-"  # not correct


def match_to_aa_bp(bp, aa_bp_ivald):
    l_hits = aa_bp_ivald[bp[0]][bp[1]]
    r_hits = aa_bp_ivald[bp[2]][bp[3]]
    l_r_hits = set(x.data for x in l_hits) & set(x.data for x in r_hits)
    for aa_bp in l_r_hits:
        # make sure it's not hitting the same side twice
        lbp, rbp = sorted(((bp[0], bp[1]), (bp[2], bp[3])))
        aa_lbp, aa_rbp = sorted(((aa_bp.lchrom, aa_bp.lpos), (aa_bp.rchrom, aa_bp.rpos)))
        if lbp[0] == aa_lbp[0] and rbp[0] == aa_rbp[0]:
            if abs(lbp[1] - aa_lbp[1]) <= bp_tol and abs(rbp[1] - aa_rbp[1]) <= bp_tol:
                #TODO:  fill this out later with the directionality checks

                return str(aa_bp)

    return 0


def parse_feature_from_feature_name(feature_name):
    sname, amp_and_ind = feature_name.rsplit("_amplicon", 1)
    anum, feat_index = amp_and_ind.split("_", 1)
    feat_type, feat_num = feat_index.rsplit("_", 1)
    return sname, anum, feat_type, feat_num


def check_if_one_feature_is_unknown(fname1, fname2):
    if fname1 == "unknown" and fname2 != "unknown":
        return 1
    elif fname2 == "unknown" and fname1 != "unknown":
        return 2

    return False


def get_location_amp(segSeqD, o_chrom, o_pos):
    hits = segSeqD[o_chrom][o_pos]
    cns = [h.data for h in hits]
    max_amp = max(cns) if cns else -1
    if max_amp > 0:
        if max_amp >= 10:
            amp_level = "_high_amp"
        elif max_amp >= 5:
            amp_level = "_low_amp"
        else:
            amp_level = "_unamplified"

    else:
        amp_level = "_nonfocal_CN"

    return amp_level


def match_to_aa_feature(bp_list, aa_feat_interval_dict, segSeqD, verbose=False):
    all_linkage_types = []
    # convert to a flattened dict that maps location to feature name
    loc_to_feat = defaultdict(IntervalTree)
    for feat_name, ivald in aa_feat_interval_dict.items():
        for chrom, ivalt in ivald.items():
            for ival in ivalt:
                loc_to_feat[chrom].addi(ival.begin, ival.end, feat_name)

    # check each identified bp
    for sv in bp_list:
        link_type = ("unknown", "unknown", "unknown")  # (left feature, right feature, connection type)
        chrom1, pos1, chrom2, pos2, dir1, dir2, total_support, split_support, disc_support, outlier_stats, _, _ = sv
        hits1 = loc_to_feat[chrom1][pos1-50:pos1+50]
        hits2 = loc_to_feat[chrom2][pos2-50:pos2+50]

        parsed_hits_1 = []
        parsed_hits_2 = []

        for ph, hl in zip([parsed_hits_1, parsed_hits_2], [hits1, hits2]):
            for h in hl:
                sname, anum, feat_type, feat_num = parse_feature_from_feature_name(h.data)
                ph.append((sname, anum, feat_type, feat_num))

        # only one end maps
        if not parsed_hits_1 and not parsed_hits_2:
            if verbose:
                print("WARNING: No feature overlaps for " + str(sv))

        elif not parsed_hits_1 or not parsed_hits_2:
            ph1_exist = False
            if parsed_hits_1:
                h_chrom, h_pos = chrom1, pos1
                o_chrom, o_pos = chrom2, pos2
                ph1_exist = True
            else:
                h_chrom, h_pos = chrom2, pos2
                o_chrom, o_pos = chrom1, pos1

            for h1 in parsed_hits_1 + parsed_hits_2:
                sname, anum, feat_type, feat_num = h1
                # Case 1: ecDNA to no feature and close by
                #if segSeqD[o_chrom][o_pos] or (h_chrom == o_chrom and abs(pos1 - pos2) < 100000):
                ft1 = feat_type if ph1_exist else "None"
                ft2 = feat_type if not ph1_exist else "None"
                amp_level = get_location_amp(segSeqD, o_chrom, o_pos)
                if h_chrom == o_chrom and abs(pos1 - pos2) < 100000:
                    link_type = (ft1, ft2, "proximal" + amp_level)

                # Case 2: ecDNA to no feature and far away same chrom
                elif h_chrom == o_chrom:
                    link_type = (ft1, ft2, "long_intrachromosomal" + amp_level)

                # Case 3: non-ecDNA to no feature and far away diff chrom
                else:
                    link_type = (ft1, ft2, "interchromosomal" + amp_level)

        else:  # if parsed_hits_1 and parsed_hits_2:
            for h1 in parsed_hits_1:
                sname1, anum1, feat_type1, feat_num1 = h1
                for h2 in parsed_hits_2:
                    sname2, anum2, feat_type2, feat_num2 = h2
                    if sname1 != sname2:  # cannot compare across different samples
                        continue

                    # Case 1: both within the same ecDNA (intra-ecDNA)
                    if (anum1, feat_type1, feat_num1) == (anum2, feat_type2, feat_num2):
                        # drill down to large and small
                        link_type = (feat_type1, feat_type2, "intrafeature")

                    # Case 2: across different ecDNAs (inter-ecDNA)
                    elif feat_type1 == feat_type2 and feat_num1 != "unknown":
                        link_type = (feat_type1, feat_type2, "interfeature_intraclass")

                    # Case 3: ecDNA to non-ecDNA & not unknown (ecDNA-to-non-ecDNA)
                    elif feat_type1 != "unknown" and feat_type2 != "unknown":
                        link_type = (feat_type1, feat_type2, "interfeature_interclass")

                    # Case 4: unknown to unknown of different or same feature
                    elif feat_type1 == "unknown" and feat_type2 == "unknown":
                        link_type = (feat_type1, feat_type2, "unknown")

                    # Case 5: ecDNA to unknown of the same AA amplicon (ecDNA-to-neighboring-amp-region)
                    elif anum1 == anum2 and check_if_one_feature_is_unknown(feat_type1, feat_type2):
                        if check_if_one_feature_is_unknown(feat_type1, feat_type2) == 1:
                            # look up graph for 1
                            amp_level = get_location_amp(segSeqD, chrom1, pos1)
                        else:
                            # look up graph for 2
                            amp_level = get_location_amp(segSeqD, chrom2, pos2)

                        link_type = (feat_type1, feat_type2, "proximal" + amp_level)

                    # Case 6: ecDNA to unknown and different AA amplicon (inter-amplification)
                    elif anum1 != anum2 and check_if_one_feature_is_unknown(feat_type1, feat_type2):
                        if check_if_one_feature_is_unknown(feat_type1, feat_type2) == 1:
                            # look up graph for 1
                            amp_level = get_location_amp(segSeqD, chrom1, pos1)
                        else:
                            # look up graph for 2
                            amp_level = get_location_amp(segSeqD, chrom2, pos2)

                        if chrom1 == chrom2:
                            link_type = (feat_type1, feat_type2, "long_intrachromosomal" + amp_level)
                        else:
                            link_type = (feat_type1, feat_type2, "interchromosomal" + amp_level)

                    else:
                        link_type = (feat_type1, feat_type2, "unknown")

        all_linkage_types.append(link_type)

    return all_linkage_types


def parse_sa_cigar(cigar_str):
    """
    Parse CIGAR string from SA tag into tuples.
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
        # ref_start, ref_end = ref_end, ref_start

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


def find_split_read_evidence(read, lcD, min_sv_size=5000, min_aln_length=20, max_overlap=10, min_mapq=5, debug=False):
    """
    Identify split read evidence for SVs from a single read.

    Args:
        read: pysam alignment object of primary alignment
        lcD: low complexity region dictionary
        min_sv_size: minimum size for intrachromosomal SVs
        min_aln_length: minimum length of clip to consider
        max_overlap: maximum allowed overlap between alignments
        debug: whether to print debug information
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

    # Calculate supplementary reference length
    # supp_ref_length = sum(l for op, l in supp_cigar if op in [0, 2, 3, 7, 8])

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

    # # if supp (aka right) has the lower genome coordinate, initiate flow through that direction
    # if is_coordinate_less_than(supp_chrom, min(supp_rstart, supp_rend), primary_chrom, min(primary_rstart, primary_rend)):

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


def create_merged_mapq_windows(mate_locations, window_size=50, extend_size=10000):
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
    Now uses proper breakpoint positions based on orientation.
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
        if dir1 =="+":
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


def find_sv_supporting_reads(bam_path, chrom, start, end, lcD, min_mapq=5, min_sv_size=5000, verbose=False):
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
            split_evidence.append(split_ev)

        # Discordant pairs
        if read.query_name not in used_reads:
            disc_ev, mate_info = find_discordant_pair_evidence(read, min_sv_size, lcD, debug=verbose)
            if disc_ev:
                used_reads.add(read.query_name)
                disc_evidence.append(disc_ev)
                mate_locations.append(mate_info)

    bamfile.close()

    fixed_mapqs = 0
    unfixed_mapqs = 0
    mate_mapq_dict = get_mate_mapqs(bam_path, mate_locations)
    for disc_ev, mate_info in zip(disc_evidence, mate_locations):
        read_name = mate_info[2]
        mate_mapq = mate_mapq_dict.get(read_name, None)  # Default to None if not found
        if mate_mapq is not None:
            fixed_mapqs+=1
        else:
            unfixed_mapqs+=1

        if disc_ev[7] is None:
            disc_ev[7] = mate_mapq
        else:
            disc_ev[8] = mate_mapq

    disc_evidence = [tuple(ev) for ev in disc_evidence]
    # Return combined evidence
    print("Retrieved {} mate mapqs for {} discordant pairs".format(fixed_mapqs, fixed_mapqs+unfixed_mapqs))
    return split_evidence + disc_evidence


def get_cluster_position(positions, evidence_types, split_weight=5.0):
    """
    Calculate refined position for a cluster using weighted evidence.

    Args:
        positions: List of positions
        evidence_types: List of corresponding evidence types
        split_weight: Weight to give split read evidence vs discordant

    Returns:
        Refined position
    """
    if not positions:
        return 0

    # Calculate weighted median
    weighted_positions = []
    for pos, ev_type in zip(positions, evidence_types):
        weight = split_weight if ev_type == 'split' else 1.0
        weighted_positions.extend([pos] * int(weight))

    # Return median of weighted positions
    weighted_positions.sort()
    mid = len(weighted_positions) // 2
    return weighted_positions[mid]


# TODO: ADD VARIANCE REPORTING
def refine_cluster_positions(cluster, cluster_distance):
    """
    Refine positions for a cluster by analyzing position distribution.
    Remove outliers and calculate most likely breakpoint positions.
    Returns refined positions and statistics about outliers removed.

    Returns:
        tuple: (pos1, pos2, filtered_cluster, outlier_stats)
    """
    positions1 = [x[1] for x in cluster]
    positions2 = [x[3] for x in cluster]
    evidence_types = [x[6] for x in cluster]

    # Sort positions and find quartiles for outlier detection
    sorted_pos1 = sorted(positions1)
    sorted_pos2 = sorted(positions2)

    # Simple IQR-based outlier detection
    def get_bounds(sorted_vals):
        q1_idx = len(sorted_vals) // 4
        q3_idx = (3 * len(sorted_vals)) // 4
        q1 = sorted_vals[max(0, q1_idx)]
        q3 = sorted_vals[min(len(sorted_vals) - 1, q3_idx)]
        iqr = max(20, q3 - q1)
        return q1 - 1.5 * iqr, q3 + 1.5 * iqr

    bounds1 = get_bounds(sorted_pos1)
    bounds2 = get_bounds(sorted_pos2)

    # Track which points were filtered and why
    outlier_stats = {
        'total_points': len(cluster),
        'outliers_pos1': 0,
        'outliers_pos2': 0,
        'outliers_both': 0,
        'retained_points': 0
    }

    # Filter out outliers and get corresponding evidence types
    filtered_cluster = []
    for i, (pos1, pos2) in enumerate(zip(positions1, positions2)):
        is_outlier1 = pos1 < bounds1[0] or pos1 > bounds1[1]
        is_outlier2 = pos2 < bounds2[0] or pos2 > bounds2[1]

        if not (is_outlier1 or is_outlier2) or evidence_types[i] == 'split':
            filtered_cluster.append(cluster[i])
        else:
            if is_outlier1 and is_outlier2:
                outlier_stats['outliers_both'] += 1
            elif is_outlier1:
                outlier_stats['outliers_pos1'] += 1
            else:
                outlier_stats['outliers_pos2'] += 1

    outlier_stats['retained_points'] = len(filtered_cluster)

    # If we filtered too many points, revert to original cluster
    if len(filtered_cluster) < len(cluster) // 2:
        filtered_cluster = cluster
        outlier_stats['reverted_to_original'] = True
    else:
        outlier_stats['reverted_to_original'] = False

    # Calculate refined positions using weighted evidence
    pos1 = get_cluster_position([x[1] for x in filtered_cluster],
                                [x[6] for x in filtered_cluster])
    pos2 = get_cluster_position([x[3] for x in filtered_cluster],
                                [x[6] for x in filtered_cluster])

    return pos1, pos2, filtered_cluster, outlier_stats


def cluster_sv_evidence(sv_evidence, cluster_distance, mean_mapq_cutoff=5, verbose=False):
    """
    Improved clustering with position-based lookahead and outlier detection.
    Uses genomic distance to determine lookahead window and tracks consumed SVs.
    """
    if not sv_evidence:
        return []

    # Normalize evidence coordinates
    normalized_evidence = []
    for ev in sv_evidence:
        chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, left_mapq, right_mapq = ev
        if is_coordinate_less_than(chrom2, pos2, chrom1, pos1):
            new_dir1 = '-' if dir2 == '+' else '+'
            new_dir2 = '-' if dir1 == '+' else '+'
            normalized_evidence.append((chrom2, pos2, chrom1, pos1, new_dir1, new_dir2, ev_type, right_mapq, left_mapq))
        else:
            normalized_evidence.append((chrom1, pos1, chrom2, pos2, dir1, dir2, ev_type, left_mapq, right_mapq))

    normalized_evidence = sort_coordinate_tuples(normalized_evidence)
    used_indices = set()  # Track which SVs have been used in clusters

    # First pass: Create initial clusters
    initial_clusters = []
    i = 0
    all_outlier_stats = []

    while i < len(normalized_evidence):
        # Skip if this SV has already been used
        if i in used_indices:
            i += 1
            continue

        current_cluster = [normalized_evidence[i]]
        used_indices.add(i)
        j = i + 1

        while j < len(normalized_evidence):
            if j in used_indices:
                j += 1
                continue

            curr_sv = normalized_evidence[j]

            # Calculate current cluster statistics
            positions1 = [x[1] for x in current_cluster]
            positions2 = [x[3] for x in current_cluster]
            median1 = sorted(positions1)[len(positions1) // 2]
            median2 = sorted(positions2)[len(positions2) // 2]

            # Check if chromosomes and orientations match
            match_chroms = (curr_sv[0] == current_cluster[0][0] and
                            curr_sv[2] == current_cluster[0][2])
            match_dirs = (curr_sv[4] == current_cluster[0][4] and
                          curr_sv[5] == current_cluster[0][5])

            # If chromosome doesn't match or we're too far away, stop looking ahead
            if curr_sv[0] != current_cluster[0][0] or abs(curr_sv[1] - median1) > 2*cluster_distance:
                break

            # If we're within distance, add to cluster
            if match_chroms and match_dirs and (abs(curr_sv[1] - median1) <= cluster_distance and
                    abs(curr_sv[3] - median2) <= cluster_distance):
                current_cluster.append(curr_sv)
                used_indices.add(j)

            j += 1

        initial_clusters.append(current_cluster)
        # Find next unused index
        while i < len(normalized_evidence) and i in used_indices:
            i += 1

    # Second pass: Refine clusters
    clustered_svs = []
    for i, cluster in enumerate(initial_clusters):
        # Refine positions and filter outliers
        refined_pos1, refined_pos2, filtered_cluster, outlier_stats = refine_cluster_positions(cluster,
                                                                                               cluster_distance)

        # Add cluster index to stats
        outlier_stats['cluster_index'] = i
        outlier_stats['original_size'] = len(cluster)
        all_outlier_stats.append(outlier_stats)

        # Count split and discordant reads
        split_count = sum(1 for x in filtered_cluster if x[6] == 'split')
        disc_count = sum(1 for x in filtered_cluster if x[6] == 'discordant')
        left_mapqs = [x[7] for x in filtered_cluster if not x[7] is None]
        right_mapqs = [x[8] for x in filtered_cluster if not x[8] is None]

        left_mean_mapq = np.mean(left_mapqs) if left_mapqs else None
        right_mean_mapq = np.mean(right_mapqs) if right_mapqs else None

        if (left_mean_mapq is not None and right_mean_mapq is not None and
                (left_mean_mapq < mean_mapq_cutoff or right_mean_mapq < mean_mapq_cutoff)):
            continue

        clustered_svs.append((
            filtered_cluster[0][0],  # chrom1
            refined_pos1,  # refined pos1
            filtered_cluster[0][2],  # chrom2
            refined_pos2,  # refined pos2
            filtered_cluster[0][4],  # dir1
            filtered_cluster[0][5],  # dir2
            len(filtered_cluster),  # total support
            split_count,  # split read support
            disc_count,  # discordant pair support
            outlier_stats,  # outlier statistics
            left_mean_mapq,
            right_mean_mapq
        ))

    return clustered_svs


def identify_connections(bam_path, feat_name, feat_ivald, lcD, aa_bp_ivald, segseqD, complete_feature_ivald,
                         output_dir="discordant_profiling_logs", verbose=False):
    """
    Check intervals in a feature for SV evidence and log results.
    Includes orientation information, split/discordant read counts, and outlier statistics.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for chrom, ivalt in feat_ivald.items():
        for ival in ivalt:
            start, end = ival.begin, ival.end
            print("Searching {}:{}-{}".format(chrom, start, end))

            logfile = os.path.join(
                output_dir,
                "{feat}_{chrom}_{start}_{end}.tsv".format(
                    feat=feat_name,
                    chrom=chrom,
                    start=start,
                    end=end
                )
            )

            candidate_svs = find_sv_supporting_reads(bam_path, chrom, start, end, lcD, verbose=verbose)
            clustered_svs = cluster_sv_evidence(candidate_svs, bp_tol, verbose=verbose)

            #match to an AA SV
            all_link_features = match_to_aa_feature(clustered_svs, complete_feature_ivald, segseqD, verbose=verbose)

            with open(logfile, 'w') as f:
                header = [
                    "sv_chrom1",
                    "sv_pos1",
                    "sv_chrom2",
                    "sv_pos2",
                    "total_read_support",
                    "split_read_support",
                    "discordant_pair_support",
                    "breakpoint_distance",
                    "is_interchromosomal",
                    "orientation_string",
                    "unfiltered_cluster_size",
                    "mean_mapq",
                    "aa_bp_evidence",
                    "connection_type",
                    "left_feat",
                    "right_feat",

                ]
                f.write('\t'.join(header) + '\n')

                for sv, feat_tup_list in zip(clustered_svs, all_link_features):
                    # Unpack the expanded sv tuple (now includes outlier_stats)
                    chrom1, pos1, chrom2, pos2, dir1, dir2, total_support, split_support, disc_support, outlier_stats, left_mean_mapq, right_mean_mapq = sv

                    is_interchrom = "1" if chrom1 != chrom2 else "0"
                    distance = str(abs(pos2 - pos1)) if chrom1 == chrom2 else "NA"
                    orientation_string = "{}{}".format(dir1, dir2)
                    in_aa = str(match_to_aa_bp((chrom1, pos1, chrom2, pos2, dir1, dir2, None), aa_bp_ivald))

                    line = [
                        chrom1,  # first SV endpoint chromosome
                        str(pos1),  # first SV endpoint position
                        chrom2,  # second SV endpoint chromosome
                        str(pos2),  # second SV endpoint position
                        str(total_support),  # total supporting reads
                        str(split_support),  # split read support
                        str(disc_support),  # discordant pair support
                        distance,  # distance between breakpoints
                        is_interchrom,  # whether SV is interchromosomal
                        orientation_string,  # combined orientation string
                        str(outlier_stats['original_size']),  # original size of cluster
                        str((round(left_mean_mapq, 1) if left_mean_mapq is not None else None,
                             round(right_mean_mapq, 1) if right_mean_mapq is not None else None)),
                        in_aa,  # this SV was also reported by AA
                        str(feat_tup_list[2]),  # connection between features
                        str(feat_tup_list[0]),  # feature of left end
                        str(feat_tup_list[1]),  # feature of right end
                    ]
                    f.write('\t'.join(line) + '\n')

            print("Processed {chrom}:{start}-{end} in feature {feat}".format(
                chrom=chrom, start=start, end=end, feat=feat_name))
            print("Found {n} potential SVs".format(n=len(clustered_svs)))
            print("Results written to {log}\n".format(log=logfile))

            # Add summary of outlier statistics
            # print("Outlier Statistics Summary:")
            # total_outliers = 0
            # total_retained = 0
            # for sv in clustered_svs:
            #     stats = sv[9]  # outlier_stats is the 10th element
            #     cluster_outliers = stats['outliers_pos1'] + stats['outliers_pos2'] + stats['outliers_both']
            #     total_outliers += cluster_outliers
            #     total_retained += stats['retained_points']
            #     print(f"Cluster {stats['cluster_index']}: "
            #           f"{cluster_outliers} outliers removed, "
            #           f"{stats['retained_points']} points retained "
            #           f"({'Reverted' if stats['reverted_to_original'] else 'Not reverted'})")
            #
            # print(f"Total: {total_outliers} outliers removed, {total_retained} points retained\n")
            # sys.exit(0)


def get_large_ref_sequences(bam, min_seq_length=1_000_000):
    """
    Get dictionary of chromosome names and lengths, filtered by minimum length.

    Args:
        bam: Open pysam AlignmentFile
        min_seq_length: Minimum sequence length to include

    Returns:
        Dictionary of chromosome names to lengths
    """
    chrom_lengths = {
        ref: length for ref, length in zip(bam.references, bam.lengths)
        if length >= min_seq_length
    }

    if not chrom_lengths:
        raise ValueError(f"No reference sequences found longer than {min_seq_length:,} bp")

    return chrom_lengths


def estimate_read_length(bam, chrom, start=20000, end=120000):
    """
    Estimate read length from a sample of reads in a region.

    Args:
        bam: Open pysam AlignmentFile
        chrom: Chromosome name to sample from
        start: Start position
        end: End position

    Returns:
        Most common read length
    """
    read_lengths = []
    for read in bam.fetch(chrom, start, end):
        if not read.is_secondary and not read.is_supplementary:
            read_lengths.append(read.query_length)
            if len(read_lengths) >= 1000:
                break

    if not read_lengths:
        raise ValueError("No valid reads found to estimate read length")

    return Counter(read_lengths).most_common(1)[0][0]


def count_valid_reads(bam, chrom, start, end, min_mapq=5):
    """
    Count valid reads in a specific region.

    Args:
        bam: Open pysam AlignmentFile
        chrom: Chromosome name
        start: Start position
        end: End position
        min_mapq: Minimum mapping quality

    Returns:
        Number of valid reads in region
    """
    return sum(1 for a in bam.fetch(chrom, start, end)
               if a.mapping_quality >= min_mapq
               and not a.is_secondary
               and not a.is_supplementary)


def estimate_bam_coverage(bam_path, num_regions=10_000, region_size=10_000, min_seq_length=1_000_000):
    """
    Estimate mean sequencing coverage of a BAM file using random sampling.
    Only considers sequences longer than min_seq_length.

    Args:
        bam_path: Path to BAM file
        num_regions: Number of random regions to sample
        region_size: Size of each region to sample in bp
        min_seq_length: Minimum sequence length to consider

    Returns:
        Estimated median coverage
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    chrom_lengths = get_large_ref_sequences(bam, min_seq_length)
    read_length = estimate_read_length(bam, next(iter(chrom_lengths)))

    coverages = []
    for _ in range(num_regions):
        chrom = random.choices(
            list(chrom_lengths.keys()),
            weights=list(chrom_lengths.values())
        )[0]

        max_pos = chrom_lengths[chrom] - region_size
        start = random.randint(0, max_pos)

        read_count = count_valid_reads(bam, chrom, start, start + region_size)
        region_coverage = (read_count * read_length) / region_size
        coverages.append(region_coverage)

    bam.close()
    return float(np.median(coverages))


def get_region_coverage_stats(bam_path, chrom, start, end, read_length=None):
    """
    Get read count and coverage for a specific region.

    Args:
        bam_path: Path to BAM file
        chrom: Chromosome name
        start: Start position
        end: End position
        read_length: Optional pre-computed read length

    Returns:
        Tuple of (read_count, coverage)
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    if read_length is None:
        read_length = estimate_read_length(bam, chrom)

    read_count = count_valid_reads(bam, chrom, start, end)
    coverage = (read_count * read_length) / (end - start)

    bam.close()
    return read_count, float(coverage)


def output_region_coverage(bam_path, feature_bed_file_list, output_file, reuse_coverage, verbose=False):
    """
    Analyze coverage of regions from specified bed files and output results to a single bed file.
    First line contains overall BAM coverage as a header comment.

    Args:
        bam_path: Path to BAM file
        feature_bed_file_list: List of bed file paths to process
        output_file: Path to output bed file
        reuse_coverage: Should the output file be reused if it already exists?
        verbose: Whether to print progress messages
    """
    # First estimate overall coverage
    if reuse_coverage and os.path.exists(output_file):
        with open(output_file, 'r') as file:
            first_line = file.readline()
            print(first_line)
            return

    estimated_coverage = estimate_bam_coverage(bam_path)
    print(f"Estimated bamfile coverage: {estimated_coverage:.2f}x")
    print(f"Storing read counts in {output_file}")

    # Open output file and write header
    with open(output_file, 'w') as fout:
        # Write coverage as header comment
        fout.write(f"#Overall_Coverage\t{estimated_coverage:.2f}\n")
        # Write column headers
        fout.write("#chr\tstart\tend\tnum_reads\n")

        # Process each bed file in the provided list
        for bed_file in feature_bed_file_list:
            if verbose:
                print(f"Processing {os.path.basename(bed_file)}")

            # Read and process each region in the bed file
            with open(bed_file) as fin:
                for line in fin:
                    if line.startswith('#'):
                        continue

                    # Parse bed line
                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])

                    # Get read count for this region
                    try:
                        read_count, _ = get_region_coverage_stats(bam_path, chrom, start, end)
                        # Write to output file
                        fout.write(f"{chrom}\t{start}\t{end}\t{read_count}\n")
                    except ValueError as e:
                        if verbose:
                            print(f"Warning: Skipping region {chrom}:{start}-{end}: {str(e)}")
                        continue


def run_profiling(bam_path, feature_bed_dir, aa_output_dir, ref, coverage_reuse, verbose):
    """
    bps is a list breakpoint objects,
    segSeqD is a defaultdict that maps to interval trees using chromosome name as a key. It is structured like
    key: chromosome, value: interval tree. segseqD['chr7'][10000:25000] would return all intervals in the segSeqD overlapping chr7:10000-25000. Can also query points.
    feat_ivald maps the name of the focal amp feature to a segSeqD of the intervals in the feature
    """

    # check if aa data repo set, construct LC datatabase
    try:
        AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + ref + "/"

    except KeyError:
        sys.stderr.write("$AA_DATA_REPO not set. Please see AA installation instructions.\n")
        sys.exit(1)

    graph_files = []
    aa_files = os.listdir(aa_output_dir)
    for file in aa_files:
        if file.endswith("_graph.txt"):
            graph_files.append(aa_output_dir + file)

    if not graph_files and not any([f.endswith("_summary.txt") for f in aa_files]):
        sys.stderr.write("No AA summary file and no graph files found in {}\n".format(aa_output_dir))
        sys.exit(1)

    # infer sample name
    sample_name_set = set()
    for graph_file in graph_files:
        sample_name_set.add(os.path.basename(graph_file.rsplit("_amplicon", 1)[0]))

    if len(sample_name_set) > 1:
        sys.stderr.write("Multiple sample names found in {}. This must contain only one sample\n".format(aa_output_dir))

    sname = sample_name_set.pop()
    print("Sample name: {}".format(sname))

    output_dir = sname + "_discordant_profiling_logs/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # retrieve bed files starting with the sample name
    feature_bed_file_list = []
    for file in os.listdir(feature_bed_dir):
        if file.startswith(sname + "_amplicon") and file.endswith("_intervals.bed"):
            feature_bed_file_list.append(feature_bed_dir + file)

    if not feature_bed_file_list:
        sys.stderr.write("No feature BED files matching sample name {} in {}\n".format(sname, feature_bed_dir))

    # Perform the coverage analysis using the collected bed files
    coverage_output = output_dir + sname + "_region_coverage.tsv"
    output_region_coverage(bam_path, feature_bed_file_list, coverage_output, coverage_reuse, verbose)

    lcD, cg5D = set_lcd(AA_DATA_REPO)

    all_aa_bps = []
    amp_to_segseqd = {}
    for gfile in graph_files:
        bps, segSeqD = parse_bpg(gfile, False)
        samp_amp = os.path.basename(gfile).rsplit("_graph.txt")[0]
        amp_to_segseqd[samp_amp] = segSeqD
        print("Got " + str(len(bps)) + " breakpoints from " + os.path.basename(gfile))
        all_aa_bps.extend(bps)

    print("")
    aa_bp_ivald = build_aa_bp_ivald(all_aa_bps)

    complete_feature_ivald = defaultdict(IntervalTree)
    for feat_bed in feature_bed_file_list:
        feat_name = os.path.basename(feat_bed).rsplit("_intervals.bed")[0]
        feat_ivald = bed_to_interval_dict(feat_bed, False)
        complete_feature_ivald[feat_name] = feat_ivald

    feat_ivald_dict = {}
    for feat_bed in feature_bed_file_list:
        feat_name = os.path.basename(feat_bed).rsplit("_intervals.bed")[0]
        amplicon_pos = feat_name.rfind("amplicon")
        next_underscore = feat_name.find("_", amplicon_pos)
        samp_amp = feat_name[:next_underscore]
        print("Pocessing " + feat_name)
        feat_ivald = bed_to_interval_dict(feat_bed, False)
        feat_ivald_dict[feat_name] = feat_ivald
        segSeqD = amp_to_segseqd[samp_amp]
        identify_connections(bam_path, feat_name, feat_ivald, lcD, aa_bp_ivald, segSeqD, complete_feature_ivald,
                             output_dir=output_dir, verbose=verbose)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify weak potential SV connections within, and exiting focal "
                                                 "amplification regions.")
    parser.add_argument("--bam", help="Path to BAM file of the sample.", type=str, required=True)
    parser.add_argument("--feature_bed_dir", help="Path to AC bed file directory.", type=str, required=True)
    parser.add_argument("--aa_results_dir", help="Path to AA results directory for one sample.", type=str,
                        required=True)
    parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, GRCh38 or mm10",
                        choices=["hg19", "GRCh37", "GRCh38", "mm10"], required=True)
    parser.add_argument("--no_reuse_coverage", action="store_true",
                        help="Do not re-compute information related to coverage in the BAM", default=False)
    parser.add_argument("--verbose", action="store_true", help="Verbose terminal output for debugging breakpoints.", default=False)
    args = parser.parse_args()

    if not args.aa_results_dir.endswith("/"):
        args.aa_results_dir += "/"

    if not args.feature_bed_dir.endswith("/"):
        args.feature_bed_dir += "/"

    run_profiling(args.bam, args.feature_bed_dir, args.aa_results_dir, args.ref, not args.no_reuse_coverage, args.verbose)



