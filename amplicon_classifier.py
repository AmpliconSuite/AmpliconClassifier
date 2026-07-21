#!/usr/bin/env python3

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"

import argparse
from collections import Counter
import contextlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import io
import inspect
import logging
from math import log
import operator
import os
import re
import shlex
import shutil
import subprocess
import sys

from intervaltree import IntervalTree, Interval

from ampclasslib.ac_annotation import *
from ampclasslib.ac_io import *
from ampclasslib.ac_util import ConfigVars, is_human_viral_hybrid, write_patched_graph
from ampclasslib.classification_records import (
    AmpliconInput, AmpliconRecord, ClassificationContext, Feature, collect_amplicon_records
)
from ampclasslib.config_params import load_config
from ampclasslib._version import __ampliconclassifier_version__
from ampclasslib.fan_ml import classify_from_edges as _fan_classify_from_edges
from ampclasslib.tid_filter import check_tid as _check_tid


def setup_logger(output_prefix, verbose=False):
    """Configure root logger to write to file and stderr."""
    log_file = "{}.log".format(output_prefix)
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger()


add_chr_tag = False
WORKER_CONTEXT = None
BFBARCHITECT_RECONSTRUCT = None
BFBARCHITECT_WRITE_GRAPH = None
BFBARCHITECT_WRITE_CYCLES = None
BFBARCHITECT_VISUALIZE = None
BFBARCHITECT_CENTROMERES = None
BFBARCHITECT_CENTROMERE_PATH = None
NO_AMP_INVALID_INTERNAL_CLASS = "No amp/Invalid"
NO_FSCNA_CLASS = "No-FSCNA"
INVALID_CLASS = "Invalid"
BFBARCHITECT_SOURCE = "BFBArchitect"
AC_SOURCE = "AC"
BFBARCHITECT_UNSUPPORTED_KWARGS_WARNED = set()


def get_bfbarchitect_version(bfbarchitect_module):
    return str(getattr(bfbarchitect_module, "__version__", None) or "unknown")


def is_non_feature_amplicon_class(amp_class):
    return amp_class in {NO_AMP_INVALID_INTERNAL_CLASS, NO_FSCNA_CLASS, INVALID_CLASS}


def amplicon_log_id(sample_name, amplicon_number):
    return "{}_{}".format(sample_name, amplicon_number)


def finalize_no_amp_invalid_class(amp_class, lc_filtered_cycle_count):
    if amp_class != NO_AMP_INVALID_INTERNAL_CLASS:
        return amp_class

    if lc_filtered_cycle_count:
        return INVALID_CLASS

    return NO_FSCNA_CLASS


def infer_fan_decomposition_class(amp_class, amp_profile, rearr_e):
    if amp_class not in {NO_AMP_INVALID_INTERNAL_CLASS, NO_FSCNA_CLASS}:
        return amp_class

    cyc_sig = amp_profile["Trivial cycle"] + amp_profile["Complex-cyclic"]
    if cyc_sig > 0:
        return "Cyclic"
    if amp_profile["Complex-non-cyclic"] > 0 or rearr_e > 1:
        return "Complex-non-cyclic"
    if amp_profile["Linear"] > 0:
        return "Linear"

    # FAN is a rearranged focal amplification mechanism; avoid pairing FAN+ with
    # a no-focal-amplification placeholder when AA decomposition was too weak.
    return "Complex-non-cyclic"

# ------------------------------------------------------------
# Methods to compute values used in classification
def get_size(cycle, segSeqD):
    return sum(segSeqD[abs(x)][2] - segSeqD[abs(x)][1] for x in cycle)


def flowWeightedCycleAmount(cycle, cn, segSeqD):
    # get length of cycle
    sc_length = get_size(cycle, segSeqD) / 1000.
    return sc_length * cn


def length_weighted_mean(values, weights):
    """Helper function to calculate weighted mean"""
    return sum(v * w for v, w in zip(values, weights)) / sum(weights)


def isRearranged(cycle, segSeqD):
    # check if it contains regions from multiple chroms
    chromList = [segSeqD[abs(ind)][0] for ind in cycle if ind != 0]
    if len(set(chromList)) > 1:
        return True

    max_del_size = 0
    for i in range(0, len(cycle) - 1):
        if cycle[i] == 0 or cycle[i + 1] == 0:
            continue
        if cycle[i] < 0 < cycle[i + 1] or cycle[i] > 0 > cycle[i + 1]:
            return True

        dist_diff = get_diff(cycle[i], cycle[i + 1], segSeqD)
        max_del_size = max(dist_diff, max_del_size)

    if max_del_size > ConfigVars.tot_min_del:
        return True

    return False


def tot_rearr_edges(graphf):
    rearr_e = 0
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                lbp, rbp = fields[1].split("->")
                _, lpd = lbp.rsplit(":")
                _, rpd = rbp.rsplit(":")

                lpos, ldir = int(lpd[:-1]), lpd[-1]
                rpos, rdir = int(rpd[:-1]), rpd[-1]
                if ldir == rdir:
                    rearr_e += 1

                elif abs(rpos - lpos) > ConfigVars.fb_dist_cut:
                    rearr_e += 1

    return rearr_e


def decompositionComplexity(graphf, cycleList, cycleCNs, segSeqD, feature_inds, exclude_inds):
    # construct intervaltree of valid regions
    hit_region_it = defaultdict(IntervalTree)
    for i in feature_inds:
        cycle = cycleList[i]
        for cv in cycle:
            if cv != 0:
                c, s, e = segSeqD[abs(cv)]
                hit_region_it[c].addi(s, e+1)

    hf_cut = 0.8
    totalGraphWeight = 0
    segs = 0
    h = "SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped".rsplit()
    h = [x.rstrip(',') for x in h]
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("SequenceEdge:"):
                h = [x.rstrip(',') for x in line.rstrip().rsplit()]

            if line.startswith("sequence"):
                fields = line.rsplit()
                c, s, e = fields[1].rsplit(":")[0], int(fields[1].rsplit(":")[1][:-1]), int(fields[2].rsplit(":")[1][:-1])+1
                if add_chr_tag and not c.startswith('chr'):
                    c = "chr" + c

                if not hit_region_it[c][s:e]:
                    continue

                cn = float(fields[3])

                fd = dict(zip(h, fields))
                size = float(fd['Size']) / 1000.

                segs += 1
                totalGraphWeight += (size * cn)
                # print(fields[1], size * cn, "individual weight")

            elif line.startswith("BreakpointEdge"):
                break

    # print(totalGraphWeight, "total graphweight")
    cycleWeights = []
    new_feat_inds = set()
    for ind, cycle in enumerate(cycleList):
        if ind not in exclude_inds:
            hits = False
            for cv in cycle:
                if cv != 0:
                    c, s, e = segSeqD[abs(cv)]
                    if hit_region_it[c][s:e]:
                        hits = True
                        break
            if hits:
                # print("DEBUG: cycle: ", cycle)
                wca = flowWeightedCycleAmount(cycle, cycleCNs[ind], segSeqD)
                if ind in feature_inds:
                    new_feat_inds.add(len(cycleWeights))

                cycleWeights.append(wca)

    cf = 0
    fe_ent = 0
    added_cf = 0
    cInd = 0
    # print("TGW", totalGraphWeight)
    if totalGraphWeight > 0:
        while cf < hf_cut and cInd < len(cycleWeights):
            if cInd in new_feat_inds:
                # print(cInd, cycleWeights[cInd], "cInd, weight")
                added_cf = cycleWeights[cInd] / float(totalGraphWeight)
                cf += added_cf
                # print(cInd, "<-cInd, weight->", added_cf, cf)
                if added_cf > 0:
                    fe_ent += (added_cf * log(added_cf))

            cInd += 1

        cf = round(cf, 5)
        # print(cf, "CF")
        rf = (1 - cf)
        if rf > 0:
            fu_ent = -1 * rf * log(rf)
        else:
            fu_ent = 0

    else:
        logger.warning("Warning: total graph weight <= 0")
        fu_ent = 0
        rf = 0

    # print("frac remain:",rf,"unexp ent", fu_ent, "exp ent", fe_ent)
    # print("DEBUG: Segs: ", segs)
    seg_ent = log(1.0 / segs) if segs > 0 else 0
    # print("DEBUG: ent, ", fu_ent - fe_ent - seg_ent)
    totalEntropy = max(0, fu_ent - fe_ent - seg_ent)
    decompEntropy = max(0, fu_ent - fe_ent)
    nsegEntropy = max(0, -1*seg_ent)
    return totalEntropy, decompEntropy, nsegEntropy


# Compute f (foldback fraction) from the edges in the AA graph alone
def compute_f_from_AA_graph(graphf, lcD):
    h = "SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped".rsplit()
    h = [x.rstrip(',') for x in h]
    fb_candidates = []  # (chrom, lo, hi, reads) — resolved for nesting after collection
    nonFbCount, maxCN, tot_over_min_cn, total_seq_len = 0, 0, 0, 0

    with open(graphf) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if line.startswith("SequenceEdge:"):
                h = [x.rstrip(',') for x in fields]

            if line.startswith("discordant"):
                lbp, rbp = fields[1].split("->")
                lchrom, lpd = lbp.rsplit(":")
                rchrom, rpd = rbp.rsplit(":")
                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                lpos, ldir = int(lpd[:-1]), lpd[-1]
                rpos, rdir = int(rpd[:-1]), rpd[-1]

                if lcD[lchrom][lpos] or lcD[rchrom][rpos]:
                    continue

                elif fields[0] == "discordant" and rchrom == lchrom and abs(rpos - lpos) <= 2000 and rdir == '-' and ldir == '+':
                    continue

                rSupp = int(fields[3])
                if ldir == rdir:
                    if lchrom == rchrom and abs(rpos - lpos) < ConfigVars.fb_dist_cut:
                        lo, hi = min(lpos, rpos), max(lpos, rpos)
                        fb_candidates.append((lchrom, lo, hi, rSupp))
                    else:
                        nonFbCount += rSupp

                else:
                    nonFbCount += rSupp

            elif line.startswith("sequence"):
                total_seq_len += int(fields[5])
                if not lcD[fields[1].rsplit(":")[0]].overlaps(int(fields[1].rsplit(":")[1][:-1]),
                                                              int(fields[2].rsplit(":")[1][:-1])):
                    ccn = float(fields[3])
                    fd = dict(zip(h, fields))
                    seglen = int(fd['Size'])
                    if seglen > 1500:
                        if ccn > maxCN:
                            maxCN = ccn
                        if ccn > ConfigVars.min_amp_cn:
                            tot_over_min_cn += seglen

    # Filter nested foldbacks: if a foldback's span [lo, hi] is completely
    # contained within another foldback's span on the same chromosome, it does
    # not represent an independent inversion event. Move its reads to the
    # denominator so the total read count is preserved for normalisation.
    fb_readcount, fbEdges = 0, 0
    for i, (ci, lo_i, hi_i, reads_i) in enumerate(fb_candidates):
        is_nested = any(
            ci == cj and lo_j <= lo_i and hi_i <= hi_j
            for j, (cj, lo_j, hi_j, _) in enumerate(fb_candidates)
            if j != i
        )
        if is_nested:
            nonFbCount += reads_i
        else:
            fb_readcount += reads_i
            fbEdges += 1

    # just return 0 if there isn't enough support
    if fbEdges < 1:
        return 0, 0, 0, maxCN, tot_over_min_cn, total_seq_len

    return fbEdges, fb_readcount, fb_readcount / max(1.0, float(fb_readcount + nonFbCount)), maxCN, tot_over_min_cn, total_seq_len


def nonbfb_cycles_are_ecdna(non_bfb_cycle_inds, cycleList, segSeqD, cycleCNs, graph_cns, amp_id=None):
    for ind in non_bfb_cycle_inds:
        cycle = cycleList[ind]
        length = get_size(cycle, segSeqD)

        if length > 100000 and cycleCNs[ind] > 5:
            if isRearranged(cycle, segSeqD):
                return True

            else:
                minSeg, maxSeg, (ca, pal), (cb, pbr) = min_max_cycle_posns(cycle, segSeqD)
                change_al = pos_lies_on_cn_change(ca, pal, graph_cns, amp_id=amp_id)
                change_br = pos_lies_on_cn_change(cb, pbr, graph_cns, amp_id=amp_id)
                return change_al and change_br

        elif args.ref == "GRCh38_viral" and is_human_viral_hybrid(args.ref, cycle, segSeqD):
            return True

    return False


def cycle_contains_foldback(cycle, segSeqD):
    # True if the cycle has a same-chromosome, opposite-orientation junction within fb_dist_cut
    # (a foldback). Mirrors the foldback test in cycles_file_bfb_props. Used to keep
    # foldback-containing cycles out of the ecDNA pool in BFB+ amplicons.
    segs = [x for x in cycle if x != 0]
    if len(segs) < 2:
        return False

    pairs = list(zip(segs, segs[1:]))
    if isCircular(cycle):
        pairs.append((segs[-1], segs[0]))

    for a, b in pairs:
        if a * b < 0 and segSeqD[abs(a)][0] == segSeqD[abs(b)][0]:
            d = get_diff(a, b, segSeqD)
            if d is not None and d < ConfigVars.fb_dist_cut:
                return True

    return False


# proportion of cycles with foldbacks
def cycles_file_bfb_props(cycleList, segSeqD, cycleCNs, invalidInds, graphf, graph_cns, amp_id=None):
    FB_breaks = 0.0
    distal_breaks = 0.0
    lin_breaks = 0.0

    bfb_weight = 0.0
    non_bfb_cycle_weight = 0.0
    tot_bfb_supp_cycles = 0

    non_bfb_cycle_inds = []
    bfb_cycle_inds = []

    for ind, ocycle in enumerate(cycleList):
        cycle = copy.copy(ocycle)
        if cycle[0] != 0:
            cycle.append(cycle[0])

        removed_zero_one_len_cycle = []
        for x in cycle:
            if segSeqD[abs(x)][2] - segSeqD[abs(x)][1] > 1:
                removed_zero_one_len_cycle.append(x)

        if not removed_zero_one_len_cycle:
            continue

        cycle = removed_zero_one_len_cycle

        # elide tiny (<1 kb) out-and-back excursions: a short segment whose cycle neighbors
        # reconnect on the same chromosome within fb_dist_cut (e.g. a small insertion sitting
        # at a foldback). Dropping it lets the underlying foldback be recognized as BFB-like.
        is_closed = len(cycle) > 1 and cycle[0] == cycle[-1]
        core = cycle[:-1] if is_closed else cycle
        if len(core) > 2:
            kept = []
            for i, x in enumerate(core):
                prev_x, next_x = core[i - 1], core[(i + 1) % len(core)]
                d = get_diff(prev_x, next_x, segSeqD)
                if segSeqD[abs(x)][2] - segSeqD[abs(x)][1] < 1000 and d is not None and d < ConfigVars.fb_dist_cut:
                    continue
                kept.append(x)

            if kept and len(kept) < len(core):
                cycle = (kept + [kept[0]]) if is_closed else kept

        hit_valid_sv = False
        isBFBelem = False
        illegalBFB = False
        for a, b in zip(cycle[:-1], cycle[1:]):
            # changes direction on same chrom
            diff = get_diff(a, b, segSeqD)
            aSize = get_size([a, ], segSeqD)
            bSize = get_size([b, ], segSeqD)
            if aSize < ConfigVars.minCycleSize and bSize < ConfigVars.minCycleSize:
                continue

            # check if front and back are connected via everted edge
            front_to_back_connection = amp_encompassed(cycle, segSeqD, graphf, add_chr_tag)
            if front_to_back_connection:
                illegalBFB = True
                hit_valid_sv = True

            else:
                if a * b < 0 and segSeqD[abs(a)][0] == segSeqD[abs(b)][0]:
                    hit_valid_sv = True
                    if diff is not None and diff < ConfigVars.fb_dist_cut:
                        isBFBelem = True
                        FB_breaks += cycleCNs[ind]

                    else:
                        distal_breaks += cycleCNs[ind]

                elif diff is None or diff > ConfigVars.tot_min_del:
                    hit_valid_sv = True
                    distal_breaks += cycleCNs[ind]

                if segSeqD[abs(a)][0] != segSeqD[abs(b)][0] and not (a == 0 or b == 0):
                    illegalBFB = True

        if illegalBFB:
            isBFBelem = False

        ocsize = get_size(ocycle, segSeqD)
        # trivial path
        if ocycle[0] == 0 and not hit_valid_sv and ocsize > 100000:
            lin_breaks += cycleCNs[ind]

        if isBFBelem:
            tot_bfb_supp_cycles += 1
            bfb_weight += cycleCNs[ind]
            bfb_cycle_inds.append(ind)

        elif ocycle[0] != 0 and ocsize > 30000:
            non_bfb_cycle_weight += cycleCNs[ind]
            non_bfb_cycle_inds.append(ind)

    hasEC = nonbfb_cycles_are_ecdna(non_bfb_cycle_inds, cycleList, segSeqD, cycleCNs, graph_cns, amp_id=amp_id)
    minBFBCyclesRequired = 2

    if set(bfb_cycle_inds).issubset(set(invalidInds)):
        return 0, 0, 0, False, non_bfb_cycle_inds, []

    if FB_breaks > 1.5 and tot_bfb_supp_cycles >= minBFBCyclesRequired:
        tot = float(FB_breaks + distal_breaks + lin_breaks)
        return FB_breaks / tot, distal_breaks / tot, bfb_weight / (non_bfb_cycle_weight + bfb_weight), hasEC, \
               non_bfb_cycle_inds, bfb_cycle_inds

    return 0, 0, 0, False, non_bfb_cycle_inds, bfb_cycle_inds


def check_max_cn(ec_cycle_inds, cycleList, segSeqD, graph_cns):
    for e_ind in ec_cycle_inds:
        for c_id in cycleList[e_ind]:
            chrom, l, r = segSeqD[abs(c_id)]
            if r - l < 1000:
                continue

            for i in graph_cns[chrom][l:r]:
                if i.data > ConfigVars.min_amp_cn:
                    return True

    return False


def get_amount_sigamp(ec_cycle_inds, cycleList, segSeqD, graph_cns):
    used_content = defaultdict(set)
    for e_ind in ec_cycle_inds:
        for c_id in cycleList[e_ind]:
            chrom, l, r = segSeqD[abs(c_id)]
            if not chrom:
                continue
            seg_t = IntervalTree([Interval(l, r+1)])
            olapping_low_cns = [x for x in graph_cns[chrom][l:r] if x.data < 4]
            for x in olapping_low_cns:
                seg_t.chop(x.begin, x.end)
            for x in seg_t:
                used_content[chrom] |= set(range(x.begin, x.end))

    total_sigamp = 0
    for chrom, useset in used_content.items():
        total_sigamp += len(useset)

    return total_sigamp


def segdup_cycle(cycle, segSeqD, graph_cns, circCyc, compCyc):
    # must be circular, simple cycle, and smallish, and not heavily rearranged
    if not circCyc or get_size(cycle, segSeqD) > ConfigVars.max_segdup_size or compCyc:
        return False

    # not excessively high CN and all CNs within 1
    cycle_seg_cns = []
    cycle_seg_lengths = []
    chrom = None
    for s in cycle:
        chrom, start, end = segSeqD[abs(s)]
        graph_segs = graph_cns[chrom][start:end]
        cycle_seg_cns.extend([x.data for x in graph_segs])
        cycle_seg_lengths.extend([x.end - x.begin for x in graph_segs])

    if not chrom or not cycle_seg_cns:
        logging.warning("Warning: found unexpected empty cycle!")
        return False

    # Check if CNs are within 0.5 of each other
    if max(cycle_seg_cns) - min(cycle_seg_cns) > 0.5:
        return False

    # Find the flanking bp coordinates
    min_coord = float('inf')
    max_coord = -1
    for s in cycle:
        _, start, end = segSeqD[abs(s)]
        if start < min_coord:
            min_coord = start
        if end > max_coord:
            max_coord = end

    # Get border CNs (250bp on each side)
    left_segs = graph_cns[chrom][min_coord-250:min_coord-1]
    right_segs = graph_cns[chrom][max_coord+1:max_coord+250]

    left_cns = [x.data for x in left_segs]
    left_lengths = [x.end - x.begin for x in left_segs]

    right_cns = [x.data for x in right_segs]
    right_lengths = [x.end - x.begin for x in right_segs]

    if not right_cns or not left_cns:
        return False

    if sum(left_lengths) < 100 or sum(right_lengths) < 100:
        return False

    # Check if border CNs are within 1 of each other
    if abs(length_weighted_mean(left_cns, left_lengths) - length_weighted_mean(right_cns, right_lengths)) > 1:
        return False

    # Calculate weighted means
    border_mean = length_weighted_mean(left_cns + right_cns, left_lengths + right_lengths)
    if border_mean > ConfigVars.sig_amp:
        return False

    # Calculate cycle weighted mean
    cycle_mean = length_weighted_mean(cycle_seg_cns, cycle_seg_lengths)

    # Check if ratio is approximately 2 (within 0.2) and not super high amplification
    ratio = cycle_mean / border_mean
    return ratio - 2 <= ConfigVars.segdup_max_extra_fraction and cycle_mean < ConfigVars.high_amp


def clusterECCycles(cycleList, cycleCNs, segSeqD, graph_cns, excludableCycleIndices=None):
    padding = 500000
    indices = [x for x in range(len(cycleList)) if cycleList[x][0] != 0 and x not in excludableCycleIndices]
    clusters = []
    seenSegs = set()
    total_EC_size = 0
    for ind in indices:
        cycle = cycleList[ind]
        csize = get_size(cycle, segSeqD)
        total_EC_size+=csize
        if cycleCNs[ind] < ConfigVars.min_flow and csize < ConfigVars.minCycleSize and not is_human_viral_hybrid(args.ref, cycle, segSeqD):
            continue

        cIndsToMerge = set()
        cycle_segs = set([segSeqD[abs(s_num)] for s_num in cycle])
        # s_set = set([segSeqD[abs(s_num)] for s_num in cycle])
        s_set = cycle_segs.difference(seenSegs)
        # s_set -= seenSegs
        if not s_set:
            continue

        for c_ind, clust_dict in enumerate(clusters):
            for s in cycle_segs:
                if clust_dict[s[0]][s[1] - padding:s[2] + padding]:
                    cIndsToMerge.add(c_ind)
                    break

        newClusters = []
        newClust = defaultdict(IntervalTree)
        for s in s_set:
            newClust[s[0]].addi(s[1], s[2] + 1, ind)

        for c_ind, currClust in enumerate(clusters):
            if c_ind in cIndsToMerge:
                for k, v in currClust.items():
                    for ival in v:
                        newClust[k].addi(ival.begin, ival.end, ival.data)

            else:
                newClusters.append(currClust)

        newClusters.append(newClust)
        clusters = newClusters
        seenSegs |= s_set

    indexClusters = []
    # extract only the cycle indices from each cluster and return
    for clust in clusters:
        currIndexSet = set()
        for k, v in clust.items():
            for ival in v:
                currIndexSet.add(ival.data)

        if get_amount_sigamp(currIndexSet, cycleList, segSeqD, graph_cns) > ConfigVars.anyCycContcut:
            indexClusters.append(currIndexSet)

    # remove those where the max CN is below threshold
    indexClusters = [x for x in indexClusters if check_max_cn(x, cycleList, segSeqD, graph_cns)]

    return indexClusters


# ------------------------------------------------------------
# Classifications

# returns True if the cycle is no-amp or invalid
def cycleIsNoAmpInvalid(cycle, cn, segSeqD, isSingleton, onlyCycle, maxCN, graph_cns, amp_id=None):
    # check if contains viral sequence
    if is_viral(args.ref, cycle, segSeqD):
        return False

    if not isSingleton and not onlyCycle:  # check if cycle contains more than one segment
        # check if it is a trivia cycle, and if so determine if the boundaries are on CN changes.
        # print("sig_amp is " + str(ConfigVars.sig_amp))
        if isCircular(cycle) and not isRearranged(cycle, segSeqD) and maxCN < ConfigVars.sig_amp:
            _, _, (ca, pal), (cb, pbr) = min_max_cycle_posns(cycle, segSeqD)
            change_al = pos_lies_on_cn_change(ca, pal, graph_cns, amp_id=amp_id)
            change_br = pos_lies_on_cn_change(cb, pbr, graph_cns, amp_id=amp_id)
            if not change_al and not change_br:
                return True

        # decomp strictness is 0.1 by default but can be changed by command-line arg
        # min_flow is 1.0 by default
        scale = min(ConfigVars.min_flow, maxCN * decomposition_strictness)

    # do something slightly stricter for singleton or only cycles since they seem to be less reliably real
    # these are simple and semi-arbitrary rules based on analysis of many samples
    elif maxCN >= ConfigVars.sig_amp:
        scale = min(3., maxCN * (1.25 * decomposition_strictness))
    else:
        # check if it shows any increase from background
        scale = 2.5


    # check if cycle flow is below threshold or max CN is below what is needed for a focal amp.
    if (cn <= scale) or (maxCN < ConfigVars.min_amp_cn):  # min_amp_cn is 4.5 by default but can be changed in config file
        return True

    length = get_size(cycle, segSeqD)

    # anything that did not already fail the copy number checks is returns true if the size is too small
    return length < ConfigVars.minCycleSize  # 5kbp for minCycleSize by default


def classifyConnections(cycleSet1, cycleSet2, clfs):
    cycleSet1, cycleSet2 = sorted([cycleSet1, cycleSet2], key=lambda x: len(x), reverse=True)
    csets = []
    resultDict = defaultdict(float)
    if not cycleSet2:
        for c1 in cycleSet1:
            csets.append(frozenset([clfs[c1], ]))

    else:
        for c1 in cycleSet1:
            for c2 in cycleSet2:
                csets.append(frozenset([clfs[c1], clfs[c2]]))

    distributed_edge_value = 1.0 / len(csets) if csets else 0
    for cset in csets:
        resultDict[cset] += distributed_edge_value

    return resultDict


# categories = ["No amp/Invalid", "Linear", "Trivial cycle", "Complex-non-cyclic", "Complex-cyclic"]
def classifyAmpliconProfile(amp_profile, rearr_e, totalCompCyclicCont, totCyclicCont, tot_over_min_cn, has_hybrid, force=False):
    if amp_profile["Virus"] == 1:
        return "Virus"

    cycSig = amp_profile["Trivial cycle"] + amp_profile["Complex-cyclic"]
    if ((cycSig > ConfigVars.cycCut or totalCompCyclicCont > ConfigVars.compCycContCut) and
            totCyclicCont > ConfigVars.anyCycContcut and tot_over_min_cn > ConfigVars.ampLenOverMinCN):
        return "Cyclic"

    elif has_hybrid and totCyclicCont > ConfigVars.anyCycContcut:
        return "Cyclic"

    elif amp_profile["Complex-non-cyclic"] + cycSig > ConfigVars.compCut:
        if rearr_e > 1 and tot_over_min_cn > ConfigVars.ampLenOverMinCN:
            return "Complex-non-cyclic"

        else:
            return "Linear"

    else:
        if max(amp_profile.values()) == 0:
            return "No amp/Invalid"

        elif amp_profile["No amp/Invalid"] > 0:
            if amp_profile["Linear"] / float(amp_profile["No amp/Invalid"]) > 0.25:
                if rearr_e >= 5:
                    return "Complex-non-cyclic"

                return "Linear"

        if force:
            del amp_profile["No amp/Invalid"]
            if cycSig > max(amp_profile.values()):
                return "Cyclic"

        maxCat = max(amp_profile.items(), key=operator.itemgetter(1))[0]
        return maxCat


def classifyBFB(fb_read_prop, fb_bwp, nonbfb_sig, bfb_cyc_ratio, maxCN, tot_over_min_cn):
    # print((fb, cyc_sig, nonbfb_sig, bfb_cyc_ratio, maxCN, tot_over_min_cn))
    if (fb_read_prop < ConfigVars.min_fb_read_prop or fb_bwp < ConfigVars.fb_break_weight_prop or
            maxCN < ConfigVars.min_amp_cn):
        return None

    # dominated by non-classical BFB cycles
    elif nonbfb_sig > ConfigVars.max_nonbfb_break_weight and bfb_cyc_ratio < ConfigVars.min_bfb_cycle_weight_ratio:
        return None

    # too small
    elif tot_over_min_cn < ConfigVars.fb_dist_cut:
        return None

    return "BFB"


def is_foldback_qc_artifact(fb_edges, fb_read_prop):
    return (
        fb_edges > ConfigVars.max_foldback_edges_qc or
        (fb_edges > 15 and fb_read_prop > 0.8)
    )


def is_dense_foldback_artifact(fb_edges, total_seq_len):
    # More than 2 foldback edges with density exceeding 1 per min_bp_per_foldback bp
    # indicates artifactual inversions rather than a true BFB structure.
    return (fb_edges > 2 and total_seq_len > 0 and
            fb_edges * ConfigVars.min_bp_per_foldback > total_seq_len)


def empty_bfbarchitect_summary():
    return {
        "min_score": "NA",
        "passing_region_count": "NA",
        "multiplicities": "NA",
        "regions": "NA",
        "whole_graph_used": "NA",
        "reverse_polarity_used": "NA",
        "passing_intervals": [],
        "passing_scores": [],
        "bfb_sources": "NA"
    }


def format_bfbarchitect_score(score):
    if score is None:
        return "NA"

    return "{:.2f}".format(score)


def tag_bfbarchitect_results(results, whole_graph_used, reverse_polarity_used):
    tagged = []
    for res in results or []:
        tagged_res = dict(res)
        tagged_res["_ac_whole_graph_used"] = bool(whole_graph_used)
        tagged_res["_ac_reverse_polarity_used"] = bool(reverse_polarity_used)
        tagged.append(tagged_res)

    return tagged


def _bfbarchitect_result_flag(res, key, fallback=False):
    value = res.get(key, fallback)
    if isinstance(value, str):
        return value == "True"

    return bool(value)


def _bfbarchitect_signature_accepts(sig_params, param_name):
    if not sig_params or param_name in sig_params:
        return True

    return any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig_params.values())


def _warn_unsupported_bfbarchitect_kwarg(param_name):
    if param_name in BFBARCHITECT_UNSUPPORTED_KWARGS_WARNED:
        return

    BFBARCHITECT_UNSUPPORTED_KWARGS_WARNED.add(param_name)
    active_logger = globals().get("logger", logging.getLogger(__name__))
    active_logger.warning(
        "Installed BFBArchitect reconstruct_bfb_from_graph() does not accept expected parameter '{}'; "
        "omitting it for compatibility. Check that AmpliconClassifier and BFBArchitect are from compatible "
        "prerelease revisions.".format(param_name)
    )


def _add_bfbarchitect_kwarg_if_supported(kwargs, sig_params, param_name, value):
    if _bfbarchitect_signature_accepts(sig_params, param_name):
        kwargs[param_name] = value
        return True

    _warn_unsupported_bfbarchitect_kwarg(param_name)
    return False


def summarize_bfbarchitect_results(results, whole_graph_used, reverse_polarity_used=False):
    all_scores = []
    passing_regions = 0
    multiplicities = []
    region_strings = []
    passing_intervals = []
    passing_scores = []
    best_score = None
    best_whole_graph_used = False
    best_reverse_polarity_used = False
    attempted_whole_graph = False
    attempted_reverse_polarity = False

    for res in results or []:
        result_whole_graph_used = _bfbarchitect_result_flag(res, "_ac_whole_graph_used", whole_graph_used)
        result_reverse_polarity_used = _bfbarchitect_result_flag(
            res, "_ac_reverse_polarity_used", reverse_polarity_used
        )
        attempted_whole_graph = attempted_whole_graph or result_whole_graph_used
        attempted_reverse_polarity = attempted_reverse_polarity or result_reverse_polarity_used

        raw_scores = res.get("scores", [])
        if raw_scores is None:
            raw_scores = []
        elif not isinstance(raw_scores, (list, tuple)):
            raw_scores = [raw_scores]

        scores = []
        for score in raw_scores:
            try:
                scores.append(float(score))
            except (TypeError, ValueError):
                continue

        region_min = min(scores) if scores else None
        if region_min is not None:
            all_scores.append(region_min)
            if best_score is None or region_min < best_score:
                best_score = region_min
                best_whole_graph_used = result_whole_graph_used
                best_reverse_polarity_used = result_reverse_polarity_used
            elif region_min == best_score:
                best_whole_graph_used = best_whole_graph_used or result_whole_graph_used
                best_reverse_polarity_used = best_reverse_polarity_used or result_reverse_polarity_used

            if region_min <= ConfigVars.bfbarchitect_max_score:
                passing_regions += 1
                region = res.get("region")
                if isinstance(region, (list, tuple)) and len(region) >= 3:
                    passing_intervals.append((region[0], int(region[1]), int(region[2])))
                passing_scores.append(region_min)

        multiplicity = res.get("multiplicity", "NA")
        multiplicities.append(str(multiplicity))

        region = res.get("region", "NA")
        if isinstance(region, (list, tuple)) and len(region) >= 3:
            region_label = "{}:{}-{}".format(region[0], region[1], region[2])
        else:
            region_label = str(region)

        region_strings.append("{}:{}:{}".format(region_label, format_bfbarchitect_score(region_min), multiplicity))

    if not all_scores and not region_strings:
        return empty_bfbarchitect_summary()

    if best_score is None:
        best_whole_graph_used = attempted_whole_graph
        best_reverse_polarity_used = attempted_reverse_polarity

    return {
        "min_score": format_bfbarchitect_score(best_score),
        "passing_region_count": str(passing_regions),
        "multiplicities": ";".join(multiplicities) if multiplicities else "NA",
        "regions": ";".join(region_strings) if region_strings else "NA",
        "whole_graph_used": str(best_whole_graph_used),
        "reverse_polarity_used": str(best_reverse_polarity_used),
        "passing_intervals": passing_intervals,
        "passing_scores": passing_scores,
        "bfb_sources": "NA"
    }


def _bfbarchitect_passing_results(results):
    passing = []
    for res in results or []:
        raw_scores = res.get("scores") or []
        if not isinstance(raw_scores, (list, tuple)):
            raw_scores = [raw_scores]
        valid = []
        for s in raw_scores:
            try:
                valid.append(float(s))
            except (TypeError, ValueError):
                pass
        if valid and min(valid) <= ConfigVars.bfbarchitect_max_score:
            passing.append(res)
    return passing


def _bfbarchitect_failing_regions(results):
    """Return region tuples for results that did not pass the score threshold."""
    failing = []
    for res in results or []:
        raw_scores = res.get("scores") or []
        if not isinstance(raw_scores, (list, tuple)):
            raw_scores = [raw_scores]
        valid = []
        for s in raw_scores:
            try:
                valid.append(float(s))
            except (TypeError, ValueError):
                pass
        if not valid or min(valid) > ConfigVars.bfbarchitect_max_score:
            region = res.get("region")
            if isinstance(region, (list, tuple)) and len(region) >= 3:
                failing.append(tuple(region))
    return failing


def bfbarchitect_whole_graph_preflight(graphf):
    sequence_edges = 0
    discordant_edges = 0
    foldback_edges = 0
    chrom_bp = defaultdict(int)

    with open(graphf) as infile:
        for line in infile:
            if line.startswith("sequence"):
                sequence_edges += 1
                fields = line.rstrip().rsplit()
                if len(fields) >= 3:
                    try:
                        lchrom, lpd = fields[1].rsplit(":", 1)
                        _rchrom, rpd = fields[2].rsplit(":", 1)
                        lpos = int(lpd[:-1])
                        rpos = int(rpd[:-1])
                        chrom_bp[lchrom] += max(0, rpos - lpos)
                    except (IndexError, ValueError):
                        pass
            elif line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                if len(fields) < 2:
                    continue

                discordant_edges += 1
                try:
                    lbp, rbp = fields[1].split("->")
                    lchrom, lpd = lbp.rsplit(":", 1)
                    rchrom, rpd = rbp.rsplit(":", 1)
                    lpos, ldir = int(lpd[:-1]), lpd[-1]
                    rpos, rdir = int(rpd[:-1]), rpd[-1]
                except (IndexError, ValueError):
                    continue
                if ldir == rdir and lchrom == rchrom and abs(rpos - lpos) <= ConfigVars.fb_dist_cut:
                    foldback_edges += 1

    foldback_fraction = foldback_edges / float(discordant_edges) if discordant_edges else 0.0
    complex_graph = (
        sequence_edges > ConfigVars.bfbarchitect_whole_graph_max_sequence_edges or
        discordant_edges > ConfigVars.bfbarchitect_whole_graph_max_discordant_edges
    )
    weak_foldback_signal = foldback_fraction < ConfigVars.bfbarchitect_whole_graph_min_foldback_fraction
    large_chrom_count = sum(
        1 for bp in chrom_bp.values() if bp >= ConfigVars.bfbarchitect_whole_graph_min_large_chrom_bp
    )
    multichromosomal_graph = large_chrom_count >= ConfigVars.bfbarchitect_whole_graph_max_large_chrom_count
    should_run = not ((complex_graph and weak_foldback_signal) or multichromosomal_graph)

    return should_run, {
        "sequence_edges": sequence_edges,
        "discordant_edges": discordant_edges,
        "foldback_edges": foldback_edges,
        "foldback_fraction": foldback_fraction,
        "large_chrom_count": large_chrom_count
    }


@contextlib.contextmanager
def suppress_bfbarchitect_chatter():
    previous_disable_level = logging.root.manager.disable
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    try:
        logging.disable(logging.INFO)
        with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
            yield
    finally:
        logging.disable(previous_disable_level)


def write_bfbarchitect_outputs(results, output_prefix, whole_graph_used):
    if not results:
        return

    if BFBARCHITECT_WRITE_GRAPH is None or BFBARCHITECT_WRITE_CYCLES is None or BFBARCHITECT_VISUALIZE is None:
        return

    output_dir = os.path.join(os.path.dirname(args.o), "bfbarchitect_outputs")
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        logger.warning("Could not create BFBArchitect output directory {}: {}".format(output_dir, str(e)))
        return

    for i, res in enumerate(results, 1):
        suffix = "whole_graph" if whole_graph_used else "region{}".format(i)
        region_prefix = os.path.join(output_dir, "{}_{}".format(output_prefix, suffix))
        graph_out = "{}_BFB_graph.txt".format(region_prefix)
        cycles_out = "{}_BFB_cycles.txt".format(region_prefix)

        try:
            with suppress_bfbarchitect_chatter():
                BFBARCHITECT_WRITE_GRAPH(graph_out, res["new_segments"], res["svs"], res["sv_info"])
                BFBARCHITECT_WRITE_CYCLES(cycles_out, res["new_segments"], res["bfb_strings"], res["scores"],
                                          res["multiplicity"])
                run_bfbarchitect_visualize(cycles_out, graph_out, "{}_BFB".format(region_prefix))
        except Exception as e:
            logger.warning("Could not write BFBArchitect outputs for {}: {}".format(region_prefix, str(e)))


def run_bfbarchitect_visualize(cycles_file, graph_file, output_prefix):
    kwargs = {
        "cycle_file": cycles_file,
        "graph_file": graph_file,
        "cnr_file": None,
        "output_prefix": output_prefix,
        "multiple": False,
    }

    try:
        sig_params = inspect.signature(BFBARCHITECT_VISUALIZE).parameters
    except (TypeError, ValueError):
        sig_params = {}

    if "centromere" in sig_params:
        kwargs["centromere"] = BFBARCHITECT_CENTROMERE_PATH
    elif "centromere_dict" in sig_params:
        kwargs["centromere_dict"] = BFBARCHITECT_CENTROMERES
    elif any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig_params.values()):
        kwargs["centromere"] = BFBARCHITECT_CENTROMERE_PATH

    return BFBARCHITECT_VISUALIZE(**kwargs)


def run_bfbarchitect_reconstruct(graphf, whole_graph, reverse_polarity=False, region=None):
    kwargs = {
        "centromere_dict": BFBARCHITECT_CENTROMERES,
        "solver": None,
        "multiple": False,
        "whole_graph": whole_graph,
        "threads": args.bfb_threads,
    }

    try:
        sig = inspect.signature(BFBARCHITECT_RECONSTRUCT)
        sig_params = sig.parameters
    except (TypeError, ValueError):
        sig_params = {}

    if region is not None:
        _add_bfbarchitect_kwarg_if_supported(kwargs, sig_params, "region", region)

    if not whole_graph:
        _add_bfbarchitect_kwarg_if_supported(
            kwargs, sig_params, "max_graph_segments", ConfigVars.bfbarchitect_max_region_segments
        )

    if ConfigVars.bfbarchitect_min_lp_bound is not None:
        _add_bfbarchitect_kwarg_if_supported(
            kwargs, sig_params, "min_lp_bound", ConfigVars.bfbarchitect_min_lp_bound
        )

    if reverse_polarity:
        _add_bfbarchitect_kwarg_if_supported(kwargs, sig_params, "reverse_polarity", True)

    def _round_scores(results):
        for res in results or []:
            raw = res.get("scores")
            if raw is None:
                res["scores"] = []
            else:
                if not isinstance(raw, (list, tuple)):
                    raw = [raw]
                res["scores"] = [round(float(s), 2) for s in raw if s is not None]
        return results

    try:
        return _round_scores(BFBARCHITECT_RECONSTRUCT(graphf, **kwargs))
    except TypeError as e:
        err = str(e)
        if "min_lp_bound" in kwargs and "min_lp_bound" in err:
            kwargs.pop("min_lp_bound")
            return _round_scores(BFBARCHITECT_RECONSTRUCT(graphf, **kwargs))
        raise


def run_bfbarchitect(graphf, fb_edges, output_prefix):
    if args.no_bfbarchitect:
        return empty_bfbarchitect_summary()

    if fb_edges == 0:
        return empty_bfbarchitect_summary()

    if BFBARCHITECT_RECONSTRUCT is None:
        return empty_bfbarchitect_summary()

    if BFBARCHITECT_CENTROMERES is None:
        return empty_bfbarchitect_summary()

    try:
        all_attempted_results = []

        # Step 1: region mode, normal polarity
        with suppress_bfbarchitect_chatter():
            results = tag_bfbarchitect_results(
                run_bfbarchitect_reconstruct(graphf, whole_graph=False, reverse_polarity=False),
                whole_graph_used=False,
                reverse_polarity_used=False
            )
        all_attempted_results.extend(results)

        passing = _bfbarchitect_passing_results(results)
        failing_regions = _bfbarchitect_failing_regions(results)

        preflight_ok = None
        graph_stats = None
        whole_graph_skipped = False

        # Step 2: if nothing passed, try whole graph with normal polarity
        if not passing:
            preflight_ok, graph_stats = bfbarchitect_whole_graph_preflight(graphf)
            if preflight_ok:
                with suppress_bfbarchitect_chatter():
                    wg_results = tag_bfbarchitect_results(
                        run_bfbarchitect_reconstruct(graphf, whole_graph=True, reverse_polarity=False),
                        whole_graph_used=True,
                        reverse_polarity_used=False
                    )
                all_attempted_results.extend(wg_results)
                wg_passing = _bfbarchitect_passing_results(wg_results)
                if wg_passing:
                    summary = summarize_bfbarchitect_results(wg_passing, whole_graph_used=True,
                                                             reverse_polarity_used=False)
                    write_bfbarchitect_outputs(wg_passing, output_prefix, whole_graph_used=True)
                    return summary
            else:
                whole_graph_skipped = True
                logger.debug(
                    "Skipping BFBArchitect whole-graph retry for {} because graph preflight failed "
                    "(sequence_edges={}, discordant_edges={}, foldback_edges={}, foldback_fraction={:.3f}, "
                    "large_chrom_count={}).".format(
                        graphf, graph_stats["sequence_edges"], graph_stats["discordant_edges"],
                        graph_stats["foldback_edges"], graph_stats["foldback_fraction"],
                        graph_stats["large_chrom_count"]
                    )
                )

        # Step 3: reverse polarity region mode for regions that failed
        rp_results = []
        if not results:
            # no candidates were found at all; let BFBArchitect re-discover with reversed polarity
            with suppress_bfbarchitect_chatter():
                rp_results = tag_bfbarchitect_results(
                    run_bfbarchitect_reconstruct(graphf, whole_graph=False, reverse_polarity=True),
                    whole_graph_used=False,
                    reverse_polarity_used=True
                )
            all_attempted_results.extend(rp_results)
        else:
            for region in failing_regions:
                with suppress_bfbarchitect_chatter():
                    rp = tag_bfbarchitect_results(
                        run_bfbarchitect_reconstruct(graphf, whole_graph=False, reverse_polarity=True,
                                                    region=region),
                        whole_graph_used=False,
                        reverse_polarity_used=True
                    )
                all_attempted_results.extend(rp)
                rp_results.extend(rp)

        rp_passing = _bfbarchitect_passing_results(rp_results)
        combined_passing = passing + rp_passing
        reverse_polarity_used = bool(rp_passing)

        # Step 4: if still nothing, try whole graph with reverse polarity
        if not combined_passing:
            if preflight_ok is None:
                preflight_ok, graph_stats = bfbarchitect_whole_graph_preflight(graphf)
            if preflight_ok:
                with suppress_bfbarchitect_chatter():
                    rp_wg_results = tag_bfbarchitect_results(
                        run_bfbarchitect_reconstruct(graphf, whole_graph=True, reverse_polarity=True),
                        whole_graph_used=True,
                        reverse_polarity_used=True
                    )
                all_attempted_results.extend(rp_wg_results)
                rp_wg_passing = _bfbarchitect_passing_results(rp_wg_results)
                if rp_wg_passing:
                    summary = summarize_bfbarchitect_results(rp_wg_passing, whole_graph_used=True,
                                                             reverse_polarity_used=True)
                    write_bfbarchitect_outputs(rp_wg_passing, output_prefix, whole_graph_used=True)
                    return summary

        display_results = combined_passing if combined_passing else all_attempted_results
        summary = summarize_bfbarchitect_results(display_results, whole_graph_used=False,
                                                  reverse_polarity_used=reverse_polarity_used)
        if whole_graph_skipped and not combined_passing:
            summary["passing_region_count"] = "0"
            summary["regions"] = "whole_graph_skipped_complex_weak_foldback"
            summary["whole_graph_used"] = "False"
        write_bfbarchitect_outputs(combined_passing, output_prefix, whole_graph_used=False)
        return summary

    except Exception as e:
        logger.warning("BFBArchitect failed for {}: {}".format(graphf, str(e)))
        return empty_bfbarchitect_summary()


def interval_dict_from_intervals(intervals):
    interval_dict = defaultdict(list)
    for chrom, start, end in intervals:
        if chrom and start < end:
            interval_dict[chrom].append((start, end))

    merge_intervals({"_": interval_dict})
    return interval_dict


def feature_intervals_overlap(left, right):
    for chrom, left_intervals in left.items():
        right_intervals = right.get(chrom, [])
        for l_start, l_end in left_intervals:
            for r_start, r_end in right_intervals:
                if max(l_start, r_start) < min(l_end, r_end):
                    return True

    return False


def add_interval_dict(target, source):
    for chrom, intervals in source.items():
        target[chrom].extend(intervals)


def graph_segment_endpoints(seq_edges):
    endpoints = defaultdict(list)
    for edge in seq_edges:
        endpoints[edge["chrom"]].append(edge["start"])
        endpoints[edge["chrom"]].append(edge["end"])

    return endpoints


def snap_position_to_nearest_endpoint(pos, endpoints):
    return min(endpoints, key=lambda x: (abs(x - pos), x))


def snap_bfbarchitect_intervals(intervals, seq_edges, add_chr_tag):
    endpoints_by_chrom = graph_segment_endpoints(seq_edges)
    snapped = []
    for chrom, start, end in intervals:
        if add_chr_tag and chrom not in endpoints_by_chrom and not chrom.startswith("chr"):
            chrom = "chr" + chrom
        endpoints = endpoints_by_chrom.get(chrom)
        if not endpoints:
            continue

        snapped_start = snap_position_to_nearest_endpoint(start, endpoints)
        snapped_end = snap_position_to_nearest_endpoint(end, endpoints)
        if snapped_start == snapped_end:
            continue

        snapped.append((chrom, min(snapped_start, snapped_end), max(snapped_start, snapped_end)))

    return snapped


def interval_overlap_bp(intervals, chrom, start, end):
    overlap = 0
    for int_start, int_end in intervals.get(chrom, []):
        overlap += max(0, min(end, int_end) - max(start, int_start))

    return overlap


def cycle_bfb_overlap_fraction(cycle, segSeqD, bfb_intervals):
    total_bp = 0
    overlap_bp = 0
    for seg in cycle:
        if seg == 0:
            continue
        chrom, start, end = segSeqD[abs(seg)]
        if not chrom or start >= end:
            continue
        total_bp += end - start
        overlap_bp += interval_overlap_bp(bfb_intervals, chrom, start, end)

    if total_bp == 0:
        return 0.0

    return overlap_bp / float(total_bp)


def cycle_indices_overlapping_bfb_intervals(cycleList, segSeqD, bfb_intervals, threshold):
    cycle_indices = set()
    for ind, cycle in enumerate(cycleList):
        if cycle_bfb_overlap_fraction(cycle, segSeqD, bfb_intervals) >= threshold:
            cycle_indices.add(ind)

    return cycle_indices


def bfb_intervals_from_cycles(cycleList, segSeqD, cycle_indices, invalidInds, graph_cns, ref):
    invalid_set = set(invalidInds)
    bfb_interval_dict = defaultdict(list)
    for b_ind in cycle_indices:
        if b_ind in invalid_set:
            continue
        for c_id in cycleList[b_ind]:
            if c_id == 0:
                continue
            chrom, start, end = segSeqD[abs(c_id)]
            if not chrom:
                continue

            seg_t = IntervalTree([Interval(start, end + 1)])
            low_cn_intervals = [
                x for x in graph_cns[chrom][start:end + 1]
                if x.data < 4 and not is_human_viral_hybrid(ref, cycleList[b_ind], segSeqD)
            ]
            for x in low_cn_intervals:
                seg_t.chop(x.begin, x.end + 1)

            for x in seg_t:
                bfb_interval_dict[chrom].append((x.begin, x.end))

    merge_intervals({"_": bfb_interval_dict})
    return bfb_interval_dict


def merge_bfb_features(features):
    merged = []
    for feature in features:
        target = None
        for existing in merged:
            if feature_intervals_overlap(existing.intervals, feature.intervals):
                target = existing
                break

        if target is None:
            merged.append(feature)
            continue

        add_interval_dict(target.intervals, feature.intervals)
        merge_intervals({"_": target.intervals})
        target.cycle_indices.update(feature.cycle_indices)
        target.sources.update(feature.sources)
        if target.score is None or (feature.score is not None and feature.score < target.score):
            target.score = feature.score
        target.metadata.setdefault("merged_features", []).append(feature.feature_id)

    for ind, feature in enumerate(merged, 1):
        feature.feature_id = "BFB_{}".format(ind)

    return merged


def build_bfb_features(bfb_cycle_inds, bfbarchitect_summary, cycleList, segSeqD, invalidInds, graphFile, graph_cns,
                       native_non_bfb_inds=None):
    # Cycles that native scoring placed in the non-BFB (ecDNA-candidate) pool must not be absorbed
    # by a BFBArchitect region just because they overlap it. Native only attributes true foldback
    # cycles to BFB; mirroring that here keeps the ecDNA call independent of which method made the
    # BFB call. BFBArchitect's reported interval is unaffected (it comes from the snapped region).
    native_non_bfb_set = set(native_non_bfb_inds) if native_non_bfb_inds else set()
    features = []
    if bfb_cycle_inds:
        native_intervals = bfb_intervals_from_cycles(
            cycleList, segSeqD, set(bfb_cycle_inds), invalidInds, graph_cns, args.ref
        )
        if native_intervals:
            features.append(Feature(
                feature_id="BFB_1",
                feature_type="BFB",
                intervals=native_intervals,
                cycle_indices=set(bfb_cycle_inds),
                sources={AC_SOURCE},
                score=None,
            ))

    passing_intervals = bfbarchitect_summary.get("passing_intervals", [])
    if passing_intervals:
        seq_edges, _sv_edges, _concordant_positions = parse_graph_edges_raw(graphFile, args.add_chr_tag)
        snapped_intervals = snap_bfbarchitect_intervals(passing_intervals, seq_edges, args.add_chr_tag)
        for ind, region in enumerate(snapped_intervals, 1):
            interval_dict = interval_dict_from_intervals([region])
            cycle_indices = cycle_indices_overlapping_bfb_intervals(
                cycleList, segSeqD, interval_dict, ConfigVars.bfbarchitect_cycle_overlap_threshold
            )
            cycle_indices -= native_non_bfb_set
            score = None
            passing_scores = bfbarchitect_summary.get("passing_scores", [])
            if ind - 1 < len(passing_scores):
                score = passing_scores[ind - 1]

            features.append(Feature(
                feature_id="BFB_{}".format(len(features) + 1),
                feature_type="BFB",
                intervals=interval_dict,
                cycle_indices=cycle_indices,
                sources={BFBARCHITECT_SOURCE},
                score=score,
                metadata={"snapped_region": region},
            ))

    return merge_bfb_features(features)


def get_bfbarchitect_centromere_path(aa_data_repo_base, ref):
    if ref == "GRCh38_viral":
        c_ref = "GRCh38"
        c_file = "GRCh38_centromere.bed"
    elif ref == "GRCh37":
        c_ref = "GRCh37"
        c_file = "human_g1k_v37_centromere.bed"
    else:
        c_ref = ref
        c_file = ref + "_centromere.bed"

    return os.path.join(aa_data_repo_base, c_ref, c_file)


def load_bfbarchitect_centromeres(aa_data_repo_base, ref):
    centromere_path = get_bfbarchitect_centromere_path(aa_data_repo_base, ref)
    if not os.path.exists(centromere_path):
        logger.warning("BFBArchitect centromere BED not found: {}. Skipping BFBArchitect.".format(centromere_path))
        return None

    centromeres = {}
    with open(centromere_path) as infile:
        for line in infile:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip().rsplit()
            if len(fields) < 3:
                continue

            chrom, start, end = fields[:3]
            centromeres[chrom] = int((int(start) + int(end)) / 2)

    if not centromeres:
        logger.warning("BFBArchitect centromere BED contained no usable intervals: {}. Skipping BFBArchitect.".format(
            centromere_path))
        return None

    return centromeres


# ------------------------------------------------------------


FEATURE_SIMILARITY_TYPES = {"ecDNA", "BFB", "Complex-non-cyclic", "Linear", "FAN"}
FILTER_SIMILARITY_TYPES = {"ecDNA", "BFB", "Complex-non-cyclic", "Linear", "FAN"}
AMPLICON_SCOPE_FILTER_TYPES = {"FAN"}
FEATURE_SIMILARITY_HEADER = [
    "Amp1", "Amp2", "SimilarityScore", "SimScorePercentile", "SimScorePvalue", "AsymmetricScore1",
    "AsymmetricScore2", "GenomicSegmentScore1", "GenomicSegmentScore2", "BreakpointScore1",
    "BreakpointScore2", "JaccardGenomicSegment", "JaccardBreakpoint", "NumSharedBPs", "Amp1NumBPs",
    "Amp2NumBPs", "AmpOverlapLen", "Amp1AmpLen", "Amp2AmpLen"
]


def feature_type_from_name(feat_name):
    return feat_name.rsplit("_", 1)[0]


def feature_number_from_name(feat_name):
    return feat_name.rsplit("_", 1)[1]


def interval_dict_to_tree(interval_dict):
    ivald = defaultdict(IntervalTree)
    for chrom, intlist in interval_dict.items():
        for start, end in intlist:
            ivald[chrom].addi(start, end)

    return ivald


def merge_interval_trees(target_ivald, source_ivald):
    for chrom, interval_tree in source_ivald.items():
        for interval in interval_tree:
            target_ivald[chrom].addi(interval.begin, interval.end)


def register_feature(feature_registry, feature_id, sample_name, amplicon_number, feat_name, graph_path, cycles_path,
                     interval_dict):
    feature_registry[feature_id] = {
        "feature_id": feature_id,
        "sample_name": sample_name,
        "amplicon_number": amplicon_number,
        "feature_type": feature_type_from_name(feat_name),
        "feature_number": feature_number_from_name(feat_name),
        "graph_path": graph_path,
        "cycles_path": cycles_path,
        "png_path": cycles_path.rsplit("_cycles.txt", 1)[0] + ".png",
        "pdf_path": cycles_path.rsplit("_cycles.txt", 1)[0] + ".pdf",
        "intervals": interval_dict,
    }


def get_cross_sample_feature_pairs(feat_to_ivald, feature_registry, amps_overlap):
    # Generator: yields overlapping cross-sample pairs lazily so the full
    # O(N^2) candidate list is never materialized for large cohorts.
    feature_items = list(feat_to_ivald.items())
    for ind1 in range(len(feature_items)):
        feature_id_1, graph_pair_1 = feature_items[ind1]
        sample_name_1 = feature_registry[feature_id_1]["sample_name"]
        for ind2 in range(ind1):
            feature_id_2, graph_pair_2 = feature_items[ind2]
            if sample_name_1 == feature_registry[feature_id_2]["sample_name"]:
                continue

            if amps_overlap(graph_pair_1[1], graph_pair_2[1]):
                yield (feature_id_1, feature_id_2)


def compute_feature_similarity_scores(feature_registry, threads=1, tmpdir=None):
    from feature_similarity import amps_overlap, compute_pairwise_feature_similarity, cn_cut

    feat_to_ivald = {}
    feat_to_graph = {}
    for feature_id, feature_info in feature_registry.items():
        if feature_info["feature_type"] not in FEATURE_SIMILARITY_TYPES:
            continue

        feat_to_ivald[feature_id] = (None, interval_dict_to_tree(feature_info["intervals"]))
        feat_to_graph[feature_id] = feature_info["graph_path"]

    candidate_pairs = get_cross_sample_feature_pairs(feat_to_ivald, feature_registry, amps_overlap)

    # Returns a single-pass generator of rows already sorted by similarity score
    # (descending), backed by an on-disk merge sort when the result set is large.
    return compute_pairwise_feature_similarity(
        feat_to_graph, feat_to_ivald, candidate_pairs, add_chr_tag=add_chr_tag, lcD=lcD, cg5D=None,
        cn_cut_value=cn_cut, min_de_value=0, threads=threads, tmpdir=tmpdir
    )


def write_feature_similarity_scores(similarity_rows, output_prefix):
    n_written = 0
    with open(output_prefix + "_feature_similarity_scores.tsv", "w") as outfile:
        outfile.write("\t".join(FEATURE_SIMILARITY_HEADER) + "\n")
        for row in similarity_rows:
            outfile.write("\t".join([str(x) for x in row]) + "\n")
            n_written += 1

    logger.info("Feature similarity score rows written: {}".format(n_written))


def read_feature_similarity_scores(output_prefix):
    # Stream the (already score-sorted) similarity TSV back as split string rows.
    # choose_similarity_filter_targets() floats the numeric fields itself, so
    # string fields are fine, and it relies only on the descending score order.
    with open(output_prefix + "_feature_similarity_scores.tsv") as infile:
        next(infile, None)  # skip header
        for line in infile:
            line = line.rstrip("\n")
            if line:
                yield line.split("\t")


def feature_filter_tuple(feature_registry, feature_id):
    feature_info = feature_registry[feature_id]
    return (
        feature_info["sample_name"],
        feature_info["amplicon_number"],
        feature_info["feature_type"],
        feature_info["feature_number"],
    )


def amplicon_filter_tuple(feature_registry, feature_id):
    feature_info = feature_registry[feature_id]
    return feature_info["sample_name"], feature_info["amplicon_number"]


def add_similarity_filter_target(feature_registry, feature_id, feats_to_filter, amps_to_filter):
    if feature_registry[feature_id]["feature_type"] in AMPLICON_SCOPE_FILTER_TYPES:
        amps_to_filter.add(amplicon_filter_tuple(feature_registry, feature_id))
    else:
        feats_to_filter.add(feature_filter_tuple(feature_registry, feature_id))


def get_amplicon_feature_ids(feature_registry, sample_name, amplicon_number):
    feature_ids = []
    for feature_id, feature_info in feature_registry.items():
        if feature_info["sample_name"] == sample_name and feature_info["amplicon_number"] == amplicon_number:
            feature_ids.append(feature_id)

    return feature_ids


def choose_similarity_filter_targets(similarity_rows, feature_registry, pval):
    feats_to_filter = set()
    amps_to_filter = set()
    filter_events = []
    for fields in similarity_rows:
        feature_id_1, feature_id_2 = fields[:2]
        if feature_id_1 not in feature_registry or feature_id_2 not in feature_registry:
            continue

        feat1 = feature_registry[feature_id_1]["feature_type"]
        feat2 = feature_registry[feature_id_2]["feature_type"]
        if feat1 not in FILTER_SIMILARITY_TYPES or feat2 not in FILTER_SIMILARITY_TYPES:
            continue

        # filter linear amp pairs >0.9 in each asymmetric score
        if feat1 == "Linear" and feat2 == "Linear":
            if float(fields[5]) > 0.9 and float(fields[6]) > 0.9:
                add_similarity_filter_target(feature_registry, feature_id_1, feats_to_filter, amps_to_filter)
                add_similarity_filter_target(feature_registry, feature_id_2, feats_to_filter, amps_to_filter)
                filter_events.append((feature_id_1, feature_id_2, "linear_asymmetric_score"))

        elif float(fields[4]) <= pval:  # at least one is not a linear, apply the p-value threshold
            add_similarity_filter_target(feature_registry, feature_id_1, feats_to_filter, amps_to_filter)
            add_similarity_filter_target(feature_registry, feature_id_2, feats_to_filter, amps_to_filter)
            filter_events.append((feature_id_1, feature_id_2, "similarity_pvalue"))

        else:
            if feat1 in AMPLICON_SCOPE_FILTER_TYPES or feat2 in AMPLICON_SCOPE_FILTER_TYPES:
                continue

            break

    return feats_to_filter, amps_to_filter, filter_events


def filter_similar_amplicons(n_files, pval, feature_registry, similarity_rows, classification_results):
    # adjust the p value cutoff based on number of input amplicons
    if pval is None:
        pval = 0.05/(max(1, n_files-1))
        logger.info("adjusted p-value cutoff set to 0.05/{}={}".format(str(n_files), str(pval)))

    else:
        logger.info("p-value cutoff set to {}".format(str(pval)))

    logger.info("\nSamples are assumed to be independent as --filter_similar was set.\nFiltering highly similar amplicons"
          " across independent samples...\n")
    feats_to_filter, amps_to_filter, filter_events = choose_similarity_filter_targets(
        similarity_rows, feature_registry, pval
    )
    chromo_event_count = sum(
        1 for event in filter_events
        if feature_registry[event[0]]["feature_type"] == "FAN"
        or feature_registry[event[1]]["feature_type"] == "FAN"
    )
    logger.info(
        "Similarity filtering threshold-passing pairs: {} total, {} involving FAN".format(
            len(filter_events), chromo_event_count
        )
    )
    logger.info(
        "Similarity filtering targets: {} feature-scope removals, {} amplicon-scope removals".format(
            len(feats_to_filter), len(amps_to_filter)
        )
    )
    if filter_events:
        logger.info("Similarity filtering triggering pairs:")
        for feature_id_1, feature_id_2, reason in filter_events:
            logger.info("{}\t{}\t{}".format(feature_id_1, feature_id_2, reason))

    if not feats_to_filter and not amps_to_filter:
        return

    samp_amp_to_filt_ivald = defaultdict(lambda: defaultdict(IntervalTree))
    samp_filt_set = set()
    samp_amp_filt_set = set()
    samp_amp_to_feat = defaultdict(set)
    logger.info("The following " + str(len(feats_to_filter)) + " features will be removed:")
    for x in feats_to_filter:
        full_featname = "_".join(x)
        logger.info(full_featname)
        samp_amp = x[0] + "_" + x[1]
        samp_filt_set.add(x[0])
        samp_amp_filt_set.add(samp_amp)
        merge_interval_trees(samp_amp_to_filt_ivald[samp_amp], interval_dict_to_tree(feature_registry[full_featname]["intervals"]))
        samp_amp_to_feat[samp_amp].add((x[2], x[3]))

    logger.info("The following " + str(len(amps_to_filter)) + " amplicons will be removed:")
    for sname, ampN in sorted(amps_to_filter):
        logger.info(sname + "_" + ampN)
        samp_amp = sname + "_" + ampN
        samp_filt_set.add(sname)
        samp_amp_filt_set.add(samp_amp)
        for feature_id in get_amplicon_feature_ids(feature_registry, sname, ampN):
            feature_info = feature_registry[feature_id]
            merge_interval_trees(samp_amp_to_filt_ivald[samp_amp],
                                 interval_dict_to_tree(feature_info["intervals"]))
            samp_amp_to_feat[samp_amp].add((feature_info["feature_type"], feature_info["feature_number"]))

    # now do the filtering
    '''
    The following ClassificationResults fields are updated during similarity filtering:
    ftgd_list, ftci_list, bpgi_list, fd_list, prop_list, featComplexityD,
    AMP_classifications, fan_results, and samp_to_ec_count.
    '''

    # map amplicon to the false regions

    for sname, anum, ag_dict, nc_dict in classification_results.ftgd_list:
        if sname + "_" + anum in samp_amp_filt_set:
            filter_whole_amplicon = (sname, anum) in amps_to_filter
            for feat_name in sorted(ag_dict.keys()):
                fname, fnum = feat_name.rsplit("_", 1)
                if filter_whole_amplicon or (sname, anum, fname, fnum) in feats_to_filter:
                    ag_dict[feat_name].clear()

            for feat_name in sorted(nc_dict.keys()):
                fname, fnum = feat_name.rsplit("_", 1)
                if filter_whole_amplicon or (sname, anum, fname, fnum) in feats_to_filter:
                    nc_dict[feat_name].clear()

    for ind, x in enumerate(classification_results.ftci_list):
        # annotated_cycle_outname = os.path.basename(cyclesFile).rsplit("_cycles")[0] + "_annotated_cycles.txt"
        # outname, cycleList, cycleCNs, segSeqD, bfb_cycle_inds, ecIndexClusters, invalidInds, rearrCycleInds = x
        outname = x[0]
        cycleList = x[1]
        segSeqD = x[3]
        invalidInds = x[6]
        samp_amp = outname.rsplit("_annotated_cycles.txt")[0]
        if samp_amp in samp_amp_filt_set:
            filt_ivald = samp_amp_to_filt_ivald[samp_amp]
            for cind, cyc in enumerate(cycleList):
                for x in cyc:
                    if filt_ivald[segSeqD[abs(x)][0]][segSeqD[abs(x)][1]:segSeqD[abs(x)][2]]:
                        invalidInds.append(cind)
                        break

            classification_results.ftci_list[ind][6] = invalidInds

    for ind, (sname, bpg_linelist, feature_dict, prop_dict) in enumerate(zip(
            classification_results.sampNames, classification_results.bpgi_list,
            classification_results.fd_list, classification_results.prop_list)):
        ampN = classification_results.cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
        samp_amp = sname + "_" + ampN
        if samp_amp in samp_amp_filt_set:
            filter_whole_amplicon = (sname, ampN) in amps_to_filter
            filt_ivald = samp_amp_to_filt_ivald[samp_amp]
            for bpgi_ind, bpg_line in enumerate(bpg_linelist):
                if filt_ivald[bpg_line[0]][bpg_line[1]] or filt_ivald[bpg_line[2]][bpg_line[3]]:
                    bpg_line[6] = "None"

            keys_to_del = set()
            for feat_name, curr_fd in feature_dict.items():
                fname, fnum = feat_name.rsplit("_", 1)
                if filter_whole_amplicon or (sname, ampN, fname, fnum) in feats_to_filter:
                    keys_to_del.add(feat_name)

            for k in keys_to_del:
                del feature_dict[k]
                del prop_dict[k]

    for sname, ampN, fname, fnum in feats_to_filter:
        try:
            del classification_results.featComplexityD[(sname, ampN, fname + "_" + fnum)]

        except KeyError:
            pass

    for sname, ampN in amps_to_filter:
        keys_to_del = [x for x in classification_results.featComplexityD if x[0] == sname and x[1] == ampN]
        for key in keys_to_del:
            del classification_results.featComplexityD[key]

    for ind, sname in enumerate(classification_results.sampNames):
        ampN = classification_results.cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
        samp_amp = sname + "_" + ampN
        if samp_amp in samp_amp_filt_set:
            ampClass, ecStat, bfbStat, ecAmpliconCount = classification_results.AMP_classifications[ind]
            if (sname, ampN) in amps_to_filter:
                if ecAmpliconCount:
                    classification_results.samp_to_ec_count[sname] -= ecAmpliconCount

                classification_results.AMP_classifications[ind] = (NO_FSCNA_CLASS, False, False, 0)
                if ind < len(classification_results.fan_results):
                    classification_results.fan_results[ind]["decision"] = "not_FAN"
                    classification_results.fan_results[ind]["probability"] = 0.0
                    classification_results.fan_results[ind]["filtered_by_similarity"] = True
                continue

            feats_to_remove = [x[0] for x in samp_amp_to_feat[samp_amp]]
            was_ec_or_bfb = False
            if "BFB" in feats_to_remove:
                bfbStat = False
                was_ec_or_bfb = True

            for x in feats_to_remove:
                if "ecDNA" in x:
                    ecAmpliconCount-=1
                    classification_results.samp_to_ec_count[sname]-=1

            if ecAmpliconCount == 0 and ecStat:
                ecStat = False
                was_ec_or_bfb = True

            if not ecStat and not bfbStat:
                if not was_ec_or_bfb:
                    ampClass = NO_FSCNA_CLASS

                else:
                    # TODO: Update this based on reclassification or "layered" classification.
                    ampClass = NO_FSCNA_CLASS

            classification_results.AMP_classifications[ind] = (ampClass, ecStat, bfbStat, ecAmpliconCount)


def _run_fan_from_edges(seq_edges, sv_edges, concordant_positions=None):
    try:
        return _fan_classify_from_edges(seq_edges, sv_edges, concordant_positions=concordant_positions)
    except Exception as e:
        logger.warning("FAN classifier failed: {}".format(e))
        return {"decision": "not_FAN", "probability": 0.0, "features": {}}


def _get_intervals_from_seq_edges(seq_edges):
    intervals = defaultdict(list)
    for e in seq_edges:
        if e["chrom"]:
            intervals[e["chrom"]].append((e["start"], e["end"]))
    return dict(intervals)


def find_fan_ecdna_cycles(cycleList, cycleCNs, segSeqD, graph_cns, invalidInds, bfb_cycle_inds):
    """Find circular cycles >250 kbp where every segment's graph CN exceeds the FAN ecDNA threshold."""
    exclude = set(invalidInds) | set(bfb_cycle_inds)
    qualifying = []
    for ind, cycle in enumerate(cycleList):
        if ind in exclude or cycle[0] == 0:
            continue
        if get_size(cycle, segSeqD) < 250000:
            continue
        all_above = True
        for seg_id in cycle:
            if seg_id == 0:
                continue
            chrom, start, end = segSeqD[abs(seg_id)]
            if end - start < 1000:
                continue
            seg_cns = [iv.data for iv in graph_cns[chrom][start:end]]
            if not seg_cns or max(seg_cns) <= ConfigVars.fan_ecDNA_min_cn:
                all_above = False
                break
        if all_above:
            qualifying.append(ind)
    return qualifying


def get_raw_cycle_props(cycleList, cycleCNs, segSeqD, maxCN, rearr_e, tot_over_min_cn, graph_cns, amp_id=None):
    cycleTypes = []
    cycleWeights = []
    rearrCycleInds = set()
    totalCompCyclicCont = 0
    totCyclicCont = 0
    AMP_dvaluesDict = {x: 0.0 for x in categories}
    invalidInds = []
    has_hybrid = False
    chromSet = set()
    for ind, cycle in enumerate(cycleList):
        has_hybrid = has_hybrid or is_human_viral_hybrid(args.ref, cycle, segSeqD)
        chromSet |= set([segSeqD[abs(ind)][0] for ind in cycle if ind != 0])
        is_only_decomp = (len(cycleList) == 1)
        isSingleton = True if len(cycle) == 3 and cycle[0] == 0 else False
        if cycleIsNoAmpInvalid(cycle, cycleCNs[ind], segSeqD, isSingleton, is_only_decomp, maxCN, graph_cns, amp_id=amp_id) and not args.force:
            invalidInds.append(ind)
            cycleTypes.append("No amp/Invalid")

        else:
            circCyc = isCircular(cycle)
            compCyc = isRearranged(cycle, segSeqD)
            if args.ref == "GRCh38_viral" and not any([x.startswith("chr") for x in chromSet]):
                cycleTypes.append("Virus")

            else:
                if segdup_cycle(cycle, segSeqD, graph_cns, circCyc, compCyc):
                    cycleTypes.append("No amp/Invalid")

                else:
                    if compCyc:
                        rearrCycleInds.add(ind)
                        if circCyc:
                            totalCompCyclicCont += get_size(cycle, segSeqD)

                    if circCyc:
                        totCyclicCont += get_size(cycle, segSeqD)

                    cycleTypes.append(ampDefs[(circCyc, compCyc)])

        currWt = flowWeightedCycleAmount(cycle, cycleCNs[ind], segSeqD)
        cycleWeights.append(currWt)

    totalWeight = max(sum(cycleWeights), 1)
    for i, wt in zip(cycleTypes, cycleWeights):
        AMP_dvaluesDict[i] += (wt / totalWeight)

    # anything stored in AMP_dvaluesDict prior to running classify will get used in classification
    # make sure you're not putting in other properties before here.
    ampClass = classifyAmpliconProfile(AMP_dvaluesDict, rearr_e, totalCompCyclicCont, totCyclicCont, tot_over_min_cn, has_hybrid)

    return totalCompCyclicCont, totCyclicCont, ampClass, totalWeight, AMP_dvaluesDict, invalidInds, cycleTypes, cycleWeights, rearrCycleInds


def make_amplicon_input(fpair):
    if len(fpair) <= 2:
        raise ValueError("File list not properly formatted")

    orig_sName, cyclesFile, graphFile = fpair[:3]
    sName = re.split("_amplicon[0-9]*", orig_sName)[0]
    ampN = cyclesFile.rstrip("_cycles.txt").rsplit("_")[-1]
    return AmpliconInput(orig_sName, sName, ampN, cyclesFile, graphFile)


def warn_if_chrom_names_mismatch(amplicon_input, segSeqD, context):
    if context.args.ref == "hg19" or context.args.ref == "GRCh37":
        chroms = set([x[0] for k, x in segSeqD.items() if k != 0])
        if context.args.ref == "hg19" and not any([x.startswith("chr") for x in chroms]):
            logger.warning("chrom names: " + str(chroms))
            logger.warning("Warning! --ref hg19 was set, but chromosome naming is not consistent with hg19 ('chr1', 'chr2', etc.)\n")
        elif context.args.ref == "GRCh37" and any([x.startswith("chr") for x in chroms]) and not context.add_chr_tag:
            logger.warning("chrom names: " + str(chroms))
            logger.warning("Warning! --ref GRCh37 was set, but chromosome naming is not consistent with GRCh37 ('1', '2', etc.)\n")


def process_amplicon(amplicon_input, context):
    logger.info("Classification started: {} {}".format(
        amplicon_input.sample_name, amplicon_input.amplicon_number
    ))
    segSeqD, cycleList, cycleCNs, lc_filtered_cycle_count = parseCycle(
        amplicon_input.cycles_file, amplicon_input.graph_file, context.add_chr_tag, context.lcD,
        context.patch_links
    )
    warn_if_chrom_names_mismatch(amplicon_input, segSeqD, context)
    return run_classification(
        amplicon_input, segSeqD, cycleList, cycleCNs, lc_filtered_cycle_count,
        context.lcD, context.gene_lookup, context.ncRNA_lookup
    )


def init_worker_context(context):
    global WORKER_CONTEXT, args, add_chr_tag, decomposition_strictness, patch_links
    WORKER_CONTEXT = context
    args = context.args
    add_chr_tag = context.add_chr_tag
    decomposition_strictness = context.decomposition_strictness
    patch_links = context.patch_links


def process_amplicon_worker(amplicon_input):
    if WORKER_CONTEXT is None:
        raise RuntimeError("Worker context was not initialized")

    return process_amplicon(amplicon_input, WORKER_CONTEXT)


def log_classification_progress(completed_samples, total_samples, sample_name):
    logger.info(
        "Classification progress: {}/{} samples completed (latest: {})".format(
            completed_samples, total_samples, sample_name
        )
    )


def classify_amplicons(amplicon_inputs, classification_context, jobs):
    total_amplicons = len(amplicon_inputs)
    if total_amplicons == 0:
        return []

    sample_remaining = Counter([x.sample_name for x in amplicon_inputs])
    total_samples = len(sample_remaining)
    completed_samples = 0
    completed_amplicons = 0
    amplicon_records = [None] * total_amplicons

    def record_completed(index, record):
        nonlocal completed_amplicons, completed_samples
        amplicon_input = amplicon_inputs[index]
        amplicon_records[index] = record
        completed_amplicons += 1
        sample_remaining[amplicon_input.sample_name] -= 1
        if sample_remaining[amplicon_input.sample_name] == 0:
            completed_samples += 1
            log_classification_progress(completed_samples, total_samples, amplicon_input.sample_name)

    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs, initializer=init_worker_context,
                                 initargs=(classification_context,)) as executor:
            future_to_index = {
                executor.submit(process_amplicon_worker, amplicon_input): index
                for index, amplicon_input in enumerate(amplicon_inputs)
            }
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                record_completed(index, future.result())
    else:
        for index, amplicon_input in enumerate(amplicon_inputs):
            record_completed(index, process_amplicon(amplicon_input, classification_context))

    return amplicon_records


def run_classification(amplicon_input, segSeqD, cycleList, cycleCNs, lc_filtered_cycle_count,
                       lcD, gene_lookup, ncRNA_lookup):
    sName = amplicon_input.sample_name
    ampN = amplicon_input.amplicon_number
    cyclesFile = amplicon_input.cycles_file
    graphFile = amplicon_input.graph_file

    feature_complexity = {}
    feature_registry = {}
    full_featname_to_graph = {}
    full_featname_to_intervals = {}

    graph_cns = get_graph_cns(graphFile, args.add_chr_tag)
    # first compute some properties about the foldbacks and copy numbers
    fb_edges, fb_readcount, fb_read_prop, maxCN, tot_over_min_cn, total_seq_len = compute_f_from_AA_graph(graphFile, lcD)
    rearr_e = tot_rearr_edges(graphFile)
    (totalCompCyclicCont, totCyclicCont, ampClass, totalWeight, AMP_dvaluesDict, invalidInds, cycleTypes, cycleWeights,
     rearrCycleInds) = get_raw_cycle_props(cycleList, cycleCNs, segSeqD, maxCN, rearr_e, tot_over_min_cn, graph_cns,
                                            amp_id=amplicon_log_id(sName, ampN))

    # decomposition/amplicon complexity
    totalEnt, decompEnt, nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD, range(len(cycleList)),
                                                        set())
    AMP_dvaluesDict["Amp_complexity"] = totalEnt
    AMP_dvaluesDict["Amp_decomp_complexity"] = decompEnt
    AMP_dvaluesDict["Amp_nseg_complexity"] = nEnt

    fb_bwp, nfb_bwp, bfb_cwp, bfbHasEC, non_bfb_cycle_inds, bfb_cycle_inds = cycles_file_bfb_props(cycleList, segSeqD,
        cycleCNs, invalidInds, graphFile, graph_cns, amp_id=amplicon_log_id(sName, ampN))

    # "foldback_read_prop", "BFB_bwp", "Distal_bwp", "BFB_cwp"
    AMP_dvaluesDict["foldback_read_prop"] = fb_read_prop
    AMP_dvaluesDict["BFB_bwp"] = fb_bwp
    AMP_dvaluesDict["Distal_bwp"] = nfb_bwp
    AMP_dvaluesDict["BFB_cwp"] = bfb_cwp
    bfbClass = classifyBFB(fb_read_prop, fb_bwp, nfb_bwp, bfb_cwp, maxCN, tot_over_min_cn) if fb_edges >= 2 else None

    non_fb_rearr_e = rearr_e - fb_edges
    bfb_suppressed_as_artifact = False
    qc_filter_reason = ""
    # heuristics to catch sequencing artifact samples
    if is_foldback_qc_artifact(fb_edges, fb_read_prop):
        bfb_suppressed_as_artifact = True
        qc_filter_reason = "foldback_qc_artifact"
        logger.warning(
            "{} has high foldback edge count ({}) and foldback read fraction ({:.3f}), suggesting "
            "read-orientation artifacts. Suppressing BFB calls and BFBArchitect for this amplicon. "
            "Re-run the sample with AmpliconSuite-pipeline --foldback_pair_support_min 3 or higher if "
            "artifacts persist.".format(
                amplicon_log_id(sName, ampN), fb_edges, fb_read_prop
            )
        )

    elif is_dense_foldback_artifact(fb_edges, total_seq_len):
        bfb_suppressed_as_artifact = True
        qc_filter_reason = "dense_foldback_artifact"
        logger.warning(
            "{} has {} foldback edges over {} bp of total amplicon sequence "
            "({:.2f} foldbacks/50kbp), suggesting high-density artifactual inversions. "
            "Suppressing BFB and ecDNA detection.".format(
                amplicon_log_id(sName, ampN), fb_edges, total_seq_len,
                fb_edges / (total_seq_len / ConfigVars.min_bp_per_foldback)
            )
        )

    if bfb_suppressed_as_artifact:
        bfbClass = False
        if non_fb_rearr_e >= 4 and tot_over_min_cn > ConfigVars.compCycContCut and maxCN > ConfigVars.sig_amp:
            ampClass = "Complex-non-cyclic"
        elif tot_over_min_cn > ConfigVars.compCycContCut and maxCN > ConfigVars.sig_amp:
            ampClass = "Linear"
        else:
            ampClass = INVALID_CLASS

    pre_finalize_ampClass = ampClass
    ampClass = finalize_no_amp_invalid_class(ampClass, lc_filtered_cycle_count)
    if pre_finalize_ampClass == NO_AMP_INVALID_INTERNAL_CLASS and ampClass == INVALID_CLASS:
        logger.info(
            "{} marked Invalid because {} cycle(s) were removed by low-complexity filtering.".format(
                amplicon_log_id(sName, ampN), lc_filtered_cycle_count
            )
        )

    # Parse graph edges once here; used for FAN and the post-ecStat TID check.
    _ca_seq_edges, _ca_sv_edges, _ca_concordant_positions = parse_graph_edges_raw(graphFile, args.add_chr_tag)

    if bfb_suppressed_as_artifact:
        bfbarchitect_summary = empty_bfbarchitect_summary()
    else:
        bfb_graphFile = graphFile
        if patch_links and not args.no_bfbarchitect:
            tmp_dir = os.path.join(os.path.dirname(args.o), "tmp_patched_graphs")
            os.makedirs(tmp_dir, exist_ok=True)
            tmp_gf = os.path.join(tmp_dir, os.path.basename(graphFile))
            if write_patched_graph(graphFile, patch_links, tmp_gf):
                bfb_graphFile = tmp_gf
        bfbarchitect_summary = run_bfbarchitect(bfb_graphFile, fb_edges, "{}_{}".format(sName, ampN))

    native_bfb_cycle_inds = list(bfb_cycle_inds) if bfbClass else []
    bfb_features = build_bfb_features(
        native_bfb_cycle_inds, bfbarchitect_summary, cycleList, segSeqD, invalidInds, graphFile, graph_cns,
        native_non_bfb_inds=non_bfb_cycle_inds
    )
    bfb_cycle_inds = sorted(set().union(*[feature.cycle_indices for feature in bfb_features])) if bfb_features else []
    bfb_sources = sorted(set().union(*[feature.sources for feature in bfb_features])) if bfb_features else []
    if bfbClass and AC_SOURCE not in bfb_sources:
        bfb_sources.append(AC_SOURCE)
    bfbarchitect_summary["bfb_sources"] = "|".join(sorted(bfb_sources)) if bfb_sources else "NA"
    bfb_detected = bool(bfbClass) or bool(bfb_features)
    bfbarchitect_detected = BFBARCHITECT_SOURCE in bfb_sources

    ca_result = _run_fan_from_edges(_ca_seq_edges, _ca_sv_edges, _ca_concordant_positions)
    ca_result["amplicon_intervals"] = _get_intervals_from_seq_edges(_ca_seq_edges)
    is_fan = (ca_result["decision"] == "FAN")
    if is_fan:
        ampClass = infer_fan_decomposition_class(ampClass, AMP_dvaluesDict, rearr_e)

    ecStat = False
    bfbStat = False
    qualifying_ecdna_inds = []

    if is_fan:
        # BFB detection is unchanged; ecDNA requires strict high-CN criterion
        if bfb_detected and not is_non_feature_amplicon_class(ampClass):
            bfbStat = True
        else:
            bfb_cycle_inds = []
        if not is_non_feature_amplicon_class(ampClass):
            qualifying_ecdna_inds = find_fan_ecdna_cycles(
                cycleList, cycleCNs, segSeqD, graph_cns, invalidInds, bfb_cycle_inds
            )
            if qualifying_ecdna_inds:
                ecStat = True

    else:
        if ampClass == "Cyclic" and not bfb_detected:
            ecStat = True
            bfb_cycle_inds = []

        elif bfb_detected and not is_non_feature_amplicon_class(ampClass):
            bfbStat = True
            # An ecDNA call in a BFB+ amplicon can only rest on cycles outside the BFB set
            # (bfb_cycle_inds, which BFBArchitect may have added to) that do not themselves
            # contain a foldback — a foldback-containing cycle is BFB-explained, not ecDNA.
            final_bfb_set = set(bfb_cycle_inds)
            updated_non_bfb_inds = [i for i in non_bfb_cycle_inds
                                    if i not in final_bfb_set
                                    and not cycle_contains_foldback(cycleList[i], segSeqD)]
            if nonbfb_cycles_are_ecdna(updated_non_bfb_inds, cycleList, segSeqD, cycleCNs, graph_cns,
                                       amp_id=amplicon_log_id(sName, ampN)):
                ecStat = True

        else:
            bfb_cycle_inds = []
            bfb_features = []

    # TID suppression: a tandem inverted duplication (one -- + one ++ SV, no other SVs,
    # modest inner CN) can produce a Cyclic decomposition that looks like ecDNA.
    if ecStat and not is_fan and not bfb_detected:
        tid_result = _check_tid(
            _ca_seq_edges, _ca_sv_edges, cycleList, segSeqD,
            cn_ratio_max=ConfigVars.tid_cn_ratio_max,
            cn_ratio_max_tight_fb=ConfigVars.tid_cn_ratio_max_tight_fb,
            tight_fb_size=ConfigVars.tid_tight_fb_size,
            max_foldback_span=ConfigVars.tid_max_foldback_span,
            close_endpoint_dist=ConfigVars.tid_close_endpoint_dist,
            tid_max_inner_cn=ConfigVars.tid_max_inner_cn,
        )
        if tid_result.get('pass'):
            ecStat = False
            ampClass = "Linear"
            logger.warning(
                "{} matches tandem inverted duplication (TID) criteria — suppressing ecDNA call "
                "(span={:.1f}kb, inner/bg CN ratio={:.2f}x).".format(
                    amplicon_log_id(sName, ampN),
                    tid_result['tid_span_kb'],
                    tid_result['ratio']
                )
            )

    # determine number of ecDNA present (excluding BFB cycles)
    ecIndexClusters = []
    if ecStat:
        if is_fan:
            non_qualifying = set(range(len(cycleList))) - set(qualifying_ecdna_inds)
            excludableCycleIndices = non_qualifying | set(bfb_cycle_inds) | set(invalidInds)
        else:
            excludableCycleIndices = set(bfb_cycle_inds + invalidInds)
            if bfbStat:
                # foldback-containing cycles belong to the BFB, never to a co-called ecDNA
                excludableCycleIndices |= {i for i in range(len(cycleList))
                                           if cycle_contains_foldback(cycleList[i], segSeqD)}
        ecIndexClusters = clusterECCycles(cycleList, cycleCNs, segSeqD, graph_cns, excludableCycleIndices)
        ecAmpliconCount = max(len(ecIndexClusters), 1)

    else:
        ecAmpliconCount = 0

    # if no ecDNA-like intervals were identified, update and re-call.
    if ecStat and not ecIndexClusters:
        if not bfbStat:
            remaining_classes = [NO_AMP_INVALID_INTERNAL_CLASS, "Linear", "Complex-non-cyclic"]
            remaining_scores = [AMP_dvaluesDict[x] for x in remaining_classes]
            ampClass = remaining_classes[remaining_scores.index(max(remaining_scores))]
            ampClass = finalize_no_amp_invalid_class(ampClass, lc_filtered_cycle_count)

        ecStat = False
        ecAmpliconCount = 0

    if bfbStat:
        final_bfb_sources = sorted(set().union(*[feature.sources for feature in bfb_features])) if bfb_features else []
        if bfbClass and AC_SOURCE not in final_bfb_sources:
            final_bfb_sources.append(AC_SOURCE)
        bfbarchitect_summary["bfb_sources"] = "|".join(sorted(final_bfb_sources)) if final_bfb_sources else "AC"
    else:
        bfbarchitect_summary["bfb_sources"] = "NA"

    # write entropy for each feature
    ecEntropies = []
    if ecAmpliconCount == 1 and not ecIndexClusters:
        ecEntropies.append((totalEnt, decompEnt, nEnt))

    for ecCycleList in ecIndexClusters:
        c_ex_I = bfb_cycle_inds if bfbStat else set()
        totalEnt, decompEnt, nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD, ecCycleList, c_ex_I)
        ecEntropies.append((totalEnt, decompEnt, nEnt))

    for ind, etup in enumerate(ecEntropies):
        feature_complexity[(sName, ampN, "ecDNA_" + str(ind + 1))] = etup

    if bfbStat:
        bfb_entropy_cycle_inds = bfb_cycle_inds if bfb_cycle_inds else range(len(cycleList))
        bfb_totalEnt, bfb_decompEnt, bfb_nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
                                                                        bfb_entropy_cycle_inds, set())
        feature_complexity[(sName, ampN, "BFB_1")] = (bfb_totalEnt, bfb_decompEnt, bfb_nEnt)

    if is_fan:
        fan_cycle_inds = set(range(len(cycleList))) - set(invalidInds)
        ca_totalEnt, ca_decompEnt, ca_nEnt = decompositionComplexity(
            graphFile, cycleList, cycleCNs, segSeqD, fan_cycle_inds, set()
        )
        feature_complexity[(sName, ampN, "FAN_1")] = (ca_totalEnt, ca_decompEnt, ca_nEnt)

    bpg_linelist, gseg_cn_d, other_class_c_inds, feature_dict, prop_dict, ampClass = amplicon_annotation(cycleList, segSeqD,
        bfb_cycle_inds, ecIndexClusters, invalidInds, bfbStat, ecStat, ampClass, graphFile, args.add_chr_tag, lcD,
        args.ref, ca_result.get("amplicon_intervals") if is_fan else None, bfb_features)

    trim_sname = sName.rsplit("/")[-1].rsplit("_amplicon")[0]
    for feat_name, curr_fd in feature_dict.items():
        if curr_fd:
            full_fname = trim_sname + "_" + ampN + "_" + feat_name
            full_featname_to_graph[full_fname] = graphFile
            full_featname_to_intervals[full_fname] = curr_fd
            register_feature(feature_registry, full_fname, trim_sname, ampN, feat_name, graphFile, cyclesFile, curr_fd)

    if not bfbStat and not ecStat and ampClass in {"Linear", "Complex-non-cyclic", "Virus"}:
        feature_complexity[(sName, ampN, ampClass + "_1")] = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
            other_class_c_inds, set())

    feat_to_amped_genes = get_genes_from_intervals(gene_lookup, feature_dict, gseg_cn_d)
    feat_to_amped_ncrna = get_genes_from_intervals(ncRNA_lookup, feature_dict, gseg_cn_d)
    feature_gene_entry = [sName, ampN, feat_to_amped_genes, feat_to_amped_ncrna]

    # store this additional information
    classification = (ampClass, ecStat, bfbStat, ecAmpliconCount)
    dvalues = [AMP_dvaluesDict[x] for x in categories]

    # Flag amplicons that contain viral sequence (pure Virus amplicons and human-viral
    # hybrids alike). Viral genomes are represented as non-"chr" contigs under the
    # GRCh38_viral reference; this is the same convention used by is_human_viral_hybrid.
    amp_chroms = set(x[0] for k, x in segSeqD.items() if k != 0)
    contains_viral = args.ref == "GRCh38_viral" and any(not c.startswith("chr") for c in amp_chroms)

    annotated_cycle_outname = os.path.basename(cyclesFile).rsplit("_cycles")[0] + "_annotated_cycles.txt"
    annotated_cycles_entry = [annotated_cycle_outname, cycleList, cycleCNs, segSeqD, bfb_cycle_inds, ecIndexClusters,
                              invalidInds, rearrCycleInds]

    return AmpliconRecord(
        sample_name=sName,
        amplicon_number=ampN,
        cycles_file=cyclesFile,
        graph_file=graphFile,
        feature_gene_entry=feature_gene_entry,
        annotated_cycles_entry=annotated_cycles_entry,
        breakpoint_graph_lines=bpg_linelist,
        feature_dict=feature_dict,
        prop_dict=prop_dict,
        classification=classification,
        dvalues=dvalues,
        bfbarchitect_summary=bfbarchitect_summary,
        fan_result=ca_result,
        feature_complexity=feature_complexity,
        ec_count=ecAmpliconCount,
        samp_amp_to_graph={sName + "_" + ampN: graphFile},
        feature_registry=feature_registry,
        full_featname_to_graph=full_featname_to_graph,
        full_featname_to_intervals=full_featname_to_intervals,
        qc_filter=qc_filter_reason,
        contains_viral=contains_viral,
    )


# ------------------------------------------------------------
'''
Amplicon Classes:
1) No-FSCNA
2) Linear
3) Complex-non-cyclic
4) Virus (viral episome w/out human DNA attached in amplicon)
5) BFB
6) ecDNA
7) Invalid
'''

categories = ["No amp/Invalid", "Linear", "Trivial cycle", "Complex-non-cyclic", "Complex-cyclic", "Virus",
              "foldback_read_prop", "BFB_bwp", "Distal_bwp", "BFB_cwp", "Amp_complexity", "Amp_decomp_complexity",
              "Amp_nseg_complexity"]

# ampDefs keys indicate (cyclic, complex)
ampDefs = {(False, False): "Linear", (False, True): "Complex-non-cyclic",
           (True, False): "Trivial cycle", (True, True): "Complex-cyclic"}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify AA amplicon type")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--AA_results", help="Path to location of AA output files. Can be multiple runs in a single parent directory.")
    group.add_argument("-i", "--input", help="Alternative to --AA_results if make_input.sh was already run to produce the .input file")
    group.add_argument("-c", "--cycles", help="Alternative to --AA_results for classifying a single amplicon. Path to an AA-formatted cycles file.")
    parser.add_argument("-g", "--graph", help="Path to AA-formatted graph file (required if --cycles given)")
    parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, or GRCh38.",
                        choices=["hg19", "GRCh37", "hg38", "GRCh38", "GRCh38_viral", "mm10", "GRCm38"], required=True)

    parser.add_argument("-o", help="Output filename prefix", required=True)
    parser.add_argument("--min_flow", type=float, help="Minimum flow to consider among decomposed paths (1.0).",
                        default=1.0)
    parser.add_argument("--min_size", type=float, help="Minimum cycle size (in bp) to consider as valid amplicon "
                        "(5000).")
    parser.add_argument("--decomposition_strictness", help="Value between 0 and 1 reflecting how strictly to filter "
                        "low CN decompositions (default = 0.1). Higher values filter more of the low-weight "
                        "decompositions.", type=float, default=0.1)
    parser.add_argument("--force", help="Disable No-FSCNA/Invalid class if possible", action='store_true')
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files.",
                        action='store_true')
    parser.add_argument("--report_complexity", help="[Deprecated - on by default] Compute a measure of amplicon entropy"
                        " for each amplicon.", action='store_true', default=True)
    parser.add_argument("--verbose_classification", help="Generate verbose output with raw classification scores "
                        "and extended BFBArchitect diagnostics.",
                        action='store_true')
    parser.add_argument("--no_LC_filter", help="Do not filter low-complexity cycles. Not recommended to set this flag.",
                        action='store_true', default=False)
    parser.add_argument("--exclude_bed", help="List of regions in which to ignore classification.")
    parser.add_argument("--filter_similar", help="Only use if all samples are of independent origins (not replicates "
                        "and not multi-region biopsies). Permits filtering of false-positive amps arising in multiple "
                        "independent samples based on similarity calculation", action='store_true')
    parser.add_argument("--filter_pval",
                        help="P value cutoff to use when --filter_similar is set. Default is 0.05/(n_amps-1)",
                        type=float)
    parser.add_argument("--config", help="Path to custom parameter configuration file. If not specified, "
                        "uses default config in ampclasslib directory.")
    parser.add_argument("--no_results_table", help="Skip creation of the summary results table after classification. "
                        "Results table is only created when using --AA_results or --input (not with -c/-g single amplicon mode).",
                        action='store_false', dest='make_results_table')
    parser.add_argument("--no_bfbarchitect", help="Disable BFBArchitect integration. By default, AmpliconClassifier "
                        "uses BFBArchitect when it is installed and available.", action='store_true')
    parser.add_argument("--bfb_threads", help="Number of threads for BFBArchitect ILP solver.", type=int, default=1)
    parser.add_argument("--jobs", help="Number of amplicons to classify in parallel.", type=int, default=1)
    parser.add_argument("-v", "--version", action='version', version=__ampliconclassifier_version__)
    parser.set_defaults(make_results_table=True)
    args = parser.parse_args()
    args.bfbarchitect = False

    # print("AmpliconClassifier " + __ampliconclassifier_version__)
    # print(" ".join(sys.argv))

    # Normalize and create output directory first
    if args.o.endswith("/"):
        if args.input:
            args.o += os.path.basename(args.input).rsplit(".")[0]
        else:
            args.o += "AC"

    if not args.o.startswith("/"):
        args.o = os.path.abspath(args.o)

    # Check for spaces in the final path
    if ' ' in args.o:
        sys.stderr.write("Error: Output path cannot contain spaces. Please use a path without spaces.\n")
        sys.stderr.write("Output path: {}\n".format(args.o))
        sys.exit(1)

    outdir_loc = os.path.dirname(args.o)
    if outdir_loc and not os.path.exists(outdir_loc):
        os.makedirs(outdir_loc)

    # Setup logger
    logger = setup_logger(args.o, verbose=args.verbose_classification)

    # Log version and command
    logger.info("AmpliconClassifier " + __ampliconclassifier_version__)
    logger.info("Command: %s", " ".join(sys.argv))

    # Load configuration
    try:
        config = load_config(args.config)
        logger.info("Configuration loaded from: {}".format(args.config if args.config else 'default config'))
    except Exception as e:
        logger.error("Error loading configuration: {}".format(e))
        sys.exit(1)

    set_config_vars(config, args)

    if args.make_results_table and not args.input and not args.AA_results:
        args.make_results_table = False  # results table requires multi-amplicon input mode

    if args.jobs < 1:
        logger.error("--jobs must be >= 1")
        sys.exit(1)

    if args.add_chr_tag:
        if args.ref == "GRCh38_viral":
            logger.warning("Warning: \'--add_chr_tag\' is not supported for \'GRCh38_viral\' reference.")
            args.add_chr_tag = False

        else:
            add_chr_tag = True

    # make input if an AA directory is given
    if args.AA_results:
        import ampclasslib
        src_dir = os.path.dirname(ampclasslib.__file__)
        input_file = os.path.join(outdir_loc, "AC.input")
        cmd = [os.path.join(src_dir, "make_input.sh"), args.AA_results, outdir_loc + "/AC"]
        logger.info("Generating .input file...")
        logger.info(" ".join(shlex.quote(x) for x in cmd))

        # Capture stdout and stderr
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Log stdout as info if there's any output
        if result.stdout.strip():
            for line in result.stdout.strip().split('\n'):
                logger.info(line)

        # Log stderr as error if there's any error output
        if result.stderr.strip():
            for line in result.stderr.strip().split('\n'):
                logger.error(line)

        # Check return code
        if result.returncode != 0:
            logger.error("Failed to make input file! Please ensure each graph and cycles file are present.")
            sys.exit(1)

        args.input = input_file

    if args.input:
        summary_map = os.path.splitext(args.input)[0] + "_summary_map.txt"
        if not os.path.exists(summary_map):
            logger.error("Summary map file not found with .input file. Please re-run make_input and"
                         " ensure _summary_map.txt file co-located with .input file")
            sys.exit(1)
        else:
            logger.info("Found summary map file %s", summary_map)
    else:  # one cycle + graph provided only
        bname = os.path.basename(args.cycles)
        sname = re.split("_amplicon[0-9]*", bname)[0]
        summary_map = {sname}

    if args.ref == "hg38":
        args.ref = "GRCh38"
    elif args.ref == "GRCm38":
        args.ref = "mm10"

    if not args.no_bfbarchitect and args.ref == "hg19" and args.add_chr_tag:
        logger.error("BFBArchitect cannot be used with --ref hg19 --add_chr_tag because BFBArchitect reads raw "
                     "graph chromosome names and the matching hg19/GRCh37 centromere source is ambiguous.")
        sys.exit(1)

    patch_links = read_patch_regions(args.ref)
    ncRNAFileLoc = get_ncrna_file_loc(args.ref)

    if 0 <= args.decomposition_strictness <= 1:
        decomposition_strictness = args.decomposition_strictness
    else:
        logger.error("--decomposition_strictness must be a value between 0 and 1")
        sys.exit(1)

    if args.ref == "GRCh38_viral" and args.filter_similar:
        logger.warning('--filter_similar cannot be used with --ref GRCh38_viral. disabling --filter_similar')
        args.filter_similar = False

    # check if aa data repo set, construct low-complexity (LC) datatabase
    try:
        AA_DATA_REPO_BASE = os.environ["AA_DATA_REPO"]
        AA_DATA_REPO = AA_DATA_REPO_BASE + "/" + args.ref + "/"
        fDict = {}
        with open(AA_DATA_REPO + "file_list.txt") as infile:
            for line in infile:
                fields = line.strip().rsplit()
                fDict[fields[0]] = fields[1]

        lcPath = AA_DATA_REPO + fDict["mapability_exclude_filename"]
        lcD = defaultdict(IntervalTree)
        if not args.no_LC_filter:
            lcD = buildLCDatabase(lcPath)

    except KeyError:
        logger.error("$AA_DATA_REPO not set. Please see AA installation instructions.")
        sys.exit(1)

    if not args.no_bfbarchitect:
        try:
            import bfbarchitect
            from bfbarchitect import reconstruct_bfb_from_graph, write_bfb_graph, write_bfb_cycles, visualize_BFB
            BFBARCHITECT_RECONSTRUCT = reconstruct_bfb_from_graph
            BFBARCHITECT_WRITE_GRAPH = write_bfb_graph
            BFBARCHITECT_WRITE_CYCLES = write_bfb_cycles
            BFBARCHITECT_VISUALIZE = visualize_BFB
            BFBARCHITECT_CENTROMERE_PATH = get_bfbarchitect_centromere_path(AA_DATA_REPO_BASE, args.ref)
            BFBARCHITECT_CENTROMERES = load_bfbarchitect_centromeres(AA_DATA_REPO_BASE, args.ref)
            if BFBARCHITECT_CENTROMERES:
                args.bfbarchitect = True
                logger.info("BFBArchitect {} integration enabled.".format(
                    get_bfbarchitect_version(bfbarchitect)
                ))
            else:
                logger.warning("BFBArchitect is installed but centromeres could not be loaded. Skipping BFBArchitect.")

        except ImportError:
            logger.warning("BFBArchitect is not installed. Skipping BFBArchitect integration. Use --no_bfbarchitect "
                           "to suppress this warning.")

    if args.exclude_bed:
        addtnl_lcD = buildLCDatabase(args.exclude_bed, filtsize=0)
        for k, t in addtnl_lcD.items():
            lcD[k].update(t)

    # Handle file list creation
    if not args.input:
        tempName = args.cycles.rsplit("/")[-1].rsplit(".")[0]
        flist = [[tempName, args.cycles, args.graph]]
    else:
        flist = readFlist(args.input)

    # read the gene list
    refGeneFileLoc = AA_DATA_REPO + fDict["gene_filename"]
    gene_lookup = parse_genes(refGeneFileLoc, add_chr_tag=add_chr_tag)
    ncRNA_lookup = parse_genes(ncRNAFileLoc, add_chr_tag=add_chr_tag, strip_chr_tag=(args.ref == "GRCh37"), keep_all=True)
    classification_context = ClassificationContext(
        args, add_chr_tag, decomposition_strictness, lcD, patch_links, gene_lookup, ncRNA_lookup
    )

    amplicon_inputs = []
    for i, fpair in enumerate(flist, 1):
        try:
            amplicon_input = make_amplicon_input(fpair)
        except ValueError:
            logger.info(str(fpair))
            logger.error("File list not properly formatted\n")
            sys.exit(1)

        amplicon_inputs.append(amplicon_input)

    missing_files = []
    for amplicon_input in amplicon_inputs:
        if not os.path.exists(amplicon_input.cycles_file):
            missing_files.append(amplicon_input.cycles_file)
        if not os.path.exists(amplicon_input.graph_file):
            missing_files.append(amplicon_input.graph_file)

    if args.input and isinstance(summary_map, str):
        with open(summary_map) as infile:
            for line in infile:
                fields = line.rsplit()
                if len(fields) > 1 and not os.path.exists(fields[1]):
                    missing_files.append(fields[1])

    if missing_files:
        logger.error("Found {} missing input file(s) referenced in {}:".format(
            len(missing_files), args.input if args.input else "cycles/graph arguments"))
        for f in missing_files:
            logger.error("  Missing: {}".format(f))
        sys.exit(1)

    if amplicon_inputs:
        first_amplicon = amplicon_inputs[0]
        last_amplicon = amplicon_inputs[-1]
        if len(amplicon_inputs) == 1:
            logger.info("Prepared 1 amplicon input: {} {}".format(
                first_amplicon.sample_name, first_amplicon.amplicon_number
            ))
        else:
            logger.info("Prepared {} amplicon inputs; first={} {}, last={} {}".format(
                len(amplicon_inputs),
                first_amplicon.sample_name, first_amplicon.amplicon_number,
                last_amplicon.sample_name, last_amplicon.amplicon_number
            ))

    jobs = min(args.jobs, max(1, len(amplicon_inputs)))
    if jobs > 1:
        logger.info("Classifying {} amplicons with {} worker processes".format(len(amplicon_inputs), jobs))
        if args.bfbarchitect:
            logger.info("BFBArchitect thread budget can reach {} workers * {} BFB threads = {} solver threads".format(
                jobs, args.bfb_threads, jobs * args.bfb_threads
            ))
    else:
        logger.info("Classifying {} amplicons with 1 worker process".format(len(amplicon_inputs)))

    amplicon_records = classify_amplicons(amplicon_inputs, classification_context, jobs)

    logger.info("Classification stage completed")

    qc_filtered = [(r.sample_name, r.amplicon_number, r.qc_filter)
                   for r in amplicon_records if r.qc_filter]
    if qc_filtered:
        qc_file = args.o + "_foldback_qc_filtered.txt"
        with open(qc_file, "w") as qf:
            qf.write("sample_name\tamplicon_number\tqc_filter\n")
            for sn, an, reason in qc_filtered:
                qf.write("{}\t{}\t{}\n".format(sn, an, reason))
        logger.info("Wrote {} foldback QC-filtered amplicon(s) to {}".format(
            len(qc_filtered), qc_file))

    classification_results = collect_amplicon_records(amplicon_records)

    feature_similarity_threads = jobs * max(1, args.bfb_threads) if args.bfbarchitect else jobs
    logger.info("Computing feature similarity scores...")
    feature_similarity_rows = compute_feature_similarity_scores(
        classification_results.feature_registry, threads=feature_similarity_threads,
        tmpdir=os.path.dirname(os.path.abspath(args.o))
    )
    write_feature_similarity_scores(feature_similarity_rows, args.o)
    if args.filter_similar:
        logger.info("Filtering similar amplicons...")
        # feature_similarity_rows is a single-pass generator consumed by the
        # write above; re-read the (score-sorted) TSV for the filter pass.
        filter_similar_amplicons(
            len(flist), args.filter_pval, classification_results.feature_registry,
            read_feature_similarity_scores(args.o), classification_results
        )

    logger.info("Writing outputs...")
    write_classification_results(args, classification_results, categories, summary_map)

    if args.make_results_table:
        logger.info("Creating results table...")
        try:
            import make_results_table

            classification_file = args.o + "_amplicon_classification_profiles.tsv"
            make_results_table.make_results_table(
                input_file=args.input,
                classification_file=classification_file,
                ref=args.ref
            )
            logger.info("Results table created successfully")
        except ImportError:
            logger.warning("Could not import make_results_table.py. Skipping results table creation.")
        except Exception as e:
            logger.warning("Failed to create results table: {}".format(str(e)))

    tmp_patched_dir = os.path.join(os.path.dirname(args.o), "tmp_patched_graphs")
    if os.path.isdir(tmp_patched_dir):
        shutil.rmtree(tmp_patched_dir)

    logger.info("")
    logger.info("Complete")
