#!/usr/bin/env python3

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"

import argparse
from collections import Counter, defaultdict
import os
import random
import sys

import numpy as np
import pysam
from intervaltree import IntervalTree

# Allow use as script or as library
if __name__ == "__main__":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from ac_io import parse_bpg, bed_to_interval_dict
    from ac_util import set_lcd
    from read_support_utils import *

else:
    from ampclasslib.ac_io import parse_bpg, bed_to_interval_dict
    from ampclasslib.ac_util import set_lcd
    from ampclasslib.read_support_utils import *

random.seed(0)


def build_aa_bp_ivald(aa_bp_list, window=bp_tol):
    aa_bp_ivald = defaultdict(IntervalTree)
    for bp in aa_bp_list:
        aa_bp_ivald[bp.lchrom].addi(bp.lpos-window, bp.lpos+window, bp)
        aa_bp_ivald[bp.rchrom].addi(bp.rpos - window, bp.rpos + window, bp)

    return aa_bp_ivald


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