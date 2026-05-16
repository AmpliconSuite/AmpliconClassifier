#!/usr/bin/env python3

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"

import argparse
import copy
import logging
from math import log
import operator
import os
import re
import subprocess
import sys

from intervaltree import IntervalTree, Interval

from ampclasslib.ac_annotation import *
from ampclasslib.ac_io import *
from ampclasslib.ac_util import ConfigVars
from ampclasslib.config_params import load_config
from ampclasslib.radar_plotting import *
from ampclasslib._version import __ampliconclassifier_version__
from ampclasslib.chromoauxesis_ml import classify_from_edges as _chromoauxesis_classify_from_edges


add_chr_tag = False
BFBARCHITECT_RECONSTRUCT = None
BFBARCHITECT_WRITE_GRAPH = None
BFBARCHITECT_WRITE_CYCLES = None
BFBARCHITECT_VISUALIZE = None
BFBARCHITECT_CENTROMERES = None
NO_AMP_INVALID_INTERNAL_CLASS = "No amp/Invalid"
NO_FSCNA_CLASS = "No-FSCNA"
INVALID_CLASS = "Invalid"


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
def compute_f_from_AA_graph(graphf):
    h = "SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped".rsplit()
    h = [x.rstrip(',') for x in h]
    with open(graphf) as infile:
        fb_readcount, nonFbCount, fbEdges, maxCN, tot_over_min_cn = 0, 0, 0, 0, 0
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
                        fb_readcount += rSupp
                        fbEdges += 1

                    else:
                        nonFbCount += rSupp

                else:
                    nonFbCount += rSupp

            elif line.startswith("sequence"):
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

    # just return 0 if there isn't enough support
    if fbEdges < 2:
        return 0, 0, 0, maxCN, tot_over_min_cn

    return fbEdges, fb_readcount, fb_readcount / max(1.0, float(fb_readcount + nonFbCount)), maxCN, tot_over_min_cn


def nonbfb_cycles_are_ecdna(non_bfb_cycle_inds, cycleList, segSeqD, cycleCNs, graph_cns):
    for ind in non_bfb_cycle_inds:
        cycle = cycleList[ind]
        length = get_size(cycle, segSeqD)

        if length > 100000 and cycleCNs[ind] > 5:
            if isRearranged(cycle, segSeqD):
                return True

            else:
                minSeg, maxSeg, (ca, pal), (cb, pbr) = min_max_cycle_posns(cycle, segSeqD)
                change_al = pos_lies_on_cn_change(ca, pal, graph_cns)
                change_br = pos_lies_on_cn_change(cb, pbr, graph_cns)
                return change_al and change_br

        elif args.ref == "GRCh38_viral" and is_human_viral_hybrid(args.ref, cycle, segSeqD):
            return True

    return False


# proportion of cycles with foldbacks
def cycles_file_bfb_props(cycleList, segSeqD, cycleCNs, invalidInds, graphf, graph_cns):
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

        hit_inversion = False
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

            else:
                if a * b < 0 and segSeqD[abs(a)][0] == segSeqD[abs(b)][0]:
                    hit_inversion = True
                    if diff is not None and diff < 50000:
                        isBFBelem = True
                        FB_breaks += cycleCNs[ind]

                    else:
                        distal_breaks += cycleCNs[ind]

                elif diff is None or diff > ConfigVars.tot_min_del:
                    hit_inversion = True
                    distal_breaks += cycleCNs[ind]

                if segSeqD[abs(a)][0] != segSeqD[abs(b)][0] and not (a == 0 or b == 0):
                    illegalBFB = True

        if illegalBFB:
            isBFBelem = False

        if cycle[0] == 0 and not hit_inversion and get_size(cycle, segSeqD) > 10000:
            lin_breaks += cycleCNs[ind]

        if isBFBelem:
            tot_bfb_supp_cycles += 1
            bfb_weight += cycleCNs[ind]
            bfb_cycle_inds.append(ind)

        elif cycle[0] != 0 and get_size(cycle[:-1], segSeqD) > 30000:
            non_bfb_cycle_weight += cycleCNs[ind]
            non_bfb_cycle_inds.append(ind)

    hasEC = nonbfb_cycles_are_ecdna(non_bfb_cycle_inds, cycleList, segSeqD, cycleCNs, graph_cns)
    minBFBCyclesRequired = 2

    if set(bfb_cycle_inds).issubset(set(invalidInds)):
        return 0, 0, 0, False, [], []

    if FB_breaks > 1.5 and tot_bfb_supp_cycles >= minBFBCyclesRequired:
        tot = float(FB_breaks + distal_breaks + lin_breaks)
        return FB_breaks / tot, distal_breaks / tot, bfb_weight / (non_bfb_cycle_weight + bfb_weight), hasEC, \
               non_bfb_cycle_inds, bfb_cycle_inds

    return 0, 0, 0, False, [], []


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
        sys.stderr.write("Warning: found unexpected empty cycle!\n")
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
def cycleIsNoAmpInvalid(cycle, cn, segSeqD, isSingleton, onlyCycle, maxCN, graph_cns):
    # check if contains viral sequence
    if is_viral(args.ref, cycle, segSeqD):
        return False

    if not isSingleton and not onlyCycle:  # check if cycle contains more than one segment
        # check if it is a trivia cycle, and if so determine if the boundaries are on CN changes.
        # print("sig_amp is " + str(ConfigVars.sig_amp))
        if isCircular(cycle) and not isRearranged(cycle, segSeqD) and maxCN < ConfigVars.sig_amp:
            _, _, (ca, pal), (cb, pbr) = min_max_cycle_posns(cycle, segSeqD)
            change_al = pos_lies_on_cn_change(ca, pal, graph_cns)
            change_br = pos_lies_on_cn_change(cb, pbr, graph_cns)
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


def empty_bfbarchitect_summary():
    return {
        "min_score": "NA",
        "passing_region_count": "NA",
        "multiplicities": "NA",
        "regions": "NA",
        "whole_graph_used": "NA"
    }


def format_bfbarchitect_score(score):
    if score is None:
        return "NA"

    return "{:.6g}".format(score)


def summarize_bfbarchitect_results(results, whole_graph_used):
    all_scores = []
    passing_regions = 0
    multiplicities = []
    region_strings = []

    for res in results or []:
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
            if region_min <= ConfigVars.bfbarchitect_max_score:
                passing_regions += 1

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

    return {
        "min_score": format_bfbarchitect_score(min(all_scores) if all_scores else None),
        "passing_region_count": str(passing_regions),
        "multiplicities": ";".join(multiplicities) if multiplicities else "NA",
        "regions": ";".join(region_strings) if region_strings else "NA",
        "whole_graph_used": str(whole_graph_used)
    }


def should_retry_bfbarchitect_whole_graph(results):
    if not results:
        return True

    all_scores = []
    for res in results:
        raw_scores = res.get("scores", [])
        if raw_scores is None:
            continue
        elif not isinstance(raw_scores, (list, tuple)):
            raw_scores = [raw_scores]

        for score in raw_scores:
            try:
                all_scores.append(float(score))
            except (TypeError, ValueError):
                continue

    return not all_scores or all(score > ConfigVars.bfbarchitect_max_score for score in all_scores)


def write_bfbarchitect_outputs(results, output_prefix, whole_graph_used):
    if not args.verbose_classification:
        return

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
            BFBARCHITECT_WRITE_GRAPH(graph_out, res["new_segments"], res["svs"], res["sv_info"])
            BFBARCHITECT_WRITE_CYCLES(cycles_out, res["new_segments"], res["bfb_strings"], res["scores"],
                                      res["multiplicity"])
            BFBARCHITECT_VISUALIZE(cycle_file=cycles_out, graph_file=graph_out, cnr_file=None,
                                   output_prefix="{}_BFB".format(region_prefix), multiple=False,
                                   centromere_dict=BFBARCHITECT_CENTROMERES)
        except Exception as e:
            logger.warning("Could not write BFBArchitect outputs for {}: {}".format(region_prefix, str(e)))


def run_bfbarchitect(graphf, fb_edges, output_prefix):
    if not args.bfbarchitect:
        return empty_bfbarchitect_summary()

    if fb_edges == 0:
        return empty_bfbarchitect_summary()

    if BFBARCHITECT_RECONSTRUCT is None:
        return empty_bfbarchitect_summary()

    if BFBARCHITECT_CENTROMERES is None:
        return empty_bfbarchitect_summary()

    try:
        results = BFBARCHITECT_RECONSTRUCT(graphf, centromere_dict=BFBARCHITECT_CENTROMERES, solver=None,
                                           multiple=False, whole_graph=False, threads=args.bfb_threads)
        whole_graph_used = False
        if should_retry_bfbarchitect_whole_graph(results):
            results = BFBARCHITECT_RECONSTRUCT(graphf, centromere_dict=BFBARCHITECT_CENTROMERES, solver=None,
                                               multiple=False, whole_graph=True, threads=args.bfb_threads)
            whole_graph_used = True

        summary = summarize_bfbarchitect_results(results, whole_graph_used)
        write_bfbarchitect_outputs(results, output_prefix, whole_graph_used)
        return summary

    except Exception as e:
        logger.warning("BFBArchitect failed for {}: {}".format(graphf, str(e)))
        return empty_bfbarchitect_summary()


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


def plotting():
    textCategories = ["No-FSCNA/\nInvalid", "Linear\namplification", "Trivial\ncycle", "Complex\nnon-cyclic",
                      "Complex\ncyclic", "BFB\nfoldback"]
    if args.plotstyle == "grouped":
        logger.info("plotting")
        make_classification_radar(textCategories, AMP_dvaluesList, args.o + "_amp_class", sampNames)

    elif args.plotstyle == "individual":
        logger.info("plotting")
        for a, s in zip(AMP_dvaluesList, sampNames):
            # print(textCategories, a)
            make_classification_radar(textCategories, [a[:len(textCategories)], ], args.o + "_" + s + "_amp_class",
                                      sampNames)


FEATURE_SIMILARITY_TYPES = {"ecDNA", "BFB", "Complex-non-cyclic", "Linear", "chromoauxesis"}
FILTER_SIMILARITY_TYPES = {"ecDNA", "BFB", "Complex-non-cyclic", "Linear", "chromoauxesis"}
AMPLICON_SCOPE_FILTER_TYPES = {"chromoauxesis"}
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
    pairs = []
    feature_items = list(feat_to_ivald.items())
    for ind1 in range(len(feature_items)):
        feature_id_1, graph_pair_1 = feature_items[ind1]
        sample_name_1 = feature_registry[feature_id_1]["sample_name"]
        for ind2 in range(ind1):
            feature_id_2, graph_pair_2 = feature_items[ind2]
            if sample_name_1 == feature_registry[feature_id_2]["sample_name"]:
                continue

            if amps_overlap(graph_pair_1[1], graph_pair_2[1]):
                pairs.append((feature_id_1, feature_id_2))

    return pairs


def compute_feature_similarity_scores(feature_registry):
    from feature_similarity import amps_overlap, parse_bpg, compute_similarity, cn_cut

    feat_to_ivald = {}
    for feature_id, feature_info in feature_registry.items():
        if feature_info["feature_type"] not in FEATURE_SIMILARITY_TYPES:
            continue

        feat_to_ivald[feature_id] = (None, interval_dict_to_tree(feature_info["intervals"]))

    candidate_pairs = get_cross_sample_feature_pairs(feat_to_ivald, feature_registry, amps_overlap)
    logger.info("Total of {} cross-sample pairs of features to compare for similarity.".format(len(candidate_pairs)))

    outdata = []
    for feature_id_1, feature_id_2 in candidate_pairs:
        graph_path_1 = feature_registry[feature_id_1]["graph_path"]
        graph_path_2 = feature_registry[feature_id_2]["graph_path"]
        if graph_path_1 == graph_path_2:
            continue

        graph1 = parse_bpg(graph_path_1, add_chr_tag, lcD, subset_ivald=feat_to_ivald[feature_id_1][1],
                           cn_cut=cn_cut, cg5D=None, min_de=0)
        graph2 = parse_bpg(graph_path_2, add_chr_tag, lcD, subset_ivald=feat_to_ivald[feature_id_2][1],
                           cn_cut=cn_cut, cg5D=None, min_de=0)
        if not graph1[1] or not graph2[1]:
            continue

        compute_similarity({feature_id_1: graph1, feature_id_2: graph2}, [(feature_id_1, feature_id_2)], outdata)

    outdata.sort(key=lambda x: (x[2], x[1], x[0]), reverse=True)
    return outdata


def write_feature_similarity_scores(similarity_rows, output_prefix):
    with open(output_prefix + "_feature_similarity_scores.tsv", "w") as outfile:
        outfile.write("\t".join(FEATURE_SIMILARITY_HEADER) + "\n")
        for row in similarity_rows:
            outfile.write("\t".join([str(x) for x in row]) + "\n")

    logger.info("Feature similarity score rows written: {}".format(len(similarity_rows)))


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


def filter_similar_amplicons(n_files, pval, feature_registry, similarity_rows):
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
        if feature_registry[event[0]]["feature_type"] == "chromoauxesis"
        or feature_registry[event[1]]["feature_type"] == "chromoauxesis"
    )
    logger.info(
        "Similarity filtering threshold-passing pairs: {} total, {} involving chromoauxesis".format(
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
    The following are updated during similarity filtering
    ftgd_list = []  # store list of feature gene classifications   -- removes from dictionary
    ftci_list = []  # store list of cycles file info               -- adds invalid tag to relevant cycles
    bpgi_list = []  # store list of bpg                            -- set feature of edge to "None"
    fd_list = []  # store list of feature_dicts                    -- delete feature from feature_dict
    prop_list = [] # store list of basic amplicon properties       -- delete feature from prop_dict
    featEntropyD = {}                                              -- delete feature from feature entropy dict
    AMP_classifications = []                                       -- apply amplicon de-classification logic
    samp_to_ec_count = defaultdict(int)                            -- decrement the count as needed
    '''

    # map amplicon to the false regions

    for sname, anum, ag_dict, nc_dict in ftgd_list:
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

    for ind, x in enumerate(ftci_list):
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

            ftci_list[ind][6] = invalidInds

    for ind, (sname, bpg_linelist, feature_dict, prop_dict) in enumerate(zip(sampNames, bpgi_list, fd_list, prop_list)):
        ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
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
            del featEntropyD[(sname, ampN, fname + "_" + fnum)]

        except KeyError:
            pass

    for sname, ampN in amps_to_filter:
        keys_to_del = [x for x in featEntropyD if x[0] == sname and x[1] == ampN]
        for key in keys_to_del:
            del featEntropyD[key]

    for ind, sname in enumerate(sampNames):
        ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
        samp_amp = sname + "_" + ampN
        if samp_amp in samp_amp_filt_set:
            ampClass, ecStat, bfbStat, ecAmpliconCount = AMP_classifications[ind]
            if (sname, ampN) in amps_to_filter:
                if ecAmpliconCount:
                    samp_to_ec_count[sname] -= ecAmpliconCount

                AMP_classifications[ind] = (NO_FSCNA_CLASS, False, False, 0)
                if ind < len(chromoauxesis_results):
                    chromoauxesis_results[ind]["decision"] = "not_chromoauxesis"
                    chromoauxesis_results[ind]["probability"] = 0.0
                    chromoauxesis_results[ind]["filtered_by_similarity"] = True
                continue

            feats_to_remove = [x[0] for x in samp_amp_to_feat[samp_amp]]
            was_ec_or_bfb = False
            if "BFB" in feats_to_remove:
                bfbStat = False
                was_ec_or_bfb = True

            for x in feats_to_remove:
                if "ecDNA" in x:
                    ecAmpliconCount-=1
                    samp_to_ec_count[sname]-=1

            if ecAmpliconCount == 0 and ecStat:
                ecStat = False
                was_ec_or_bfb = True

            if not ecStat and not bfbStat:
                if not was_ec_or_bfb:
                    ampClass = NO_FSCNA_CLASS

                else:
                    # TODO: Update this based on reclassification or "layered" classification.
                    ampClass = NO_FSCNA_CLASS

            AMP_classifications[ind] = (ampClass, ecStat, bfbStat, ecAmpliconCount)


def _run_chromoauxesis_from_edges(seq_edges, sv_edges):
    try:
        return _chromoauxesis_classify_from_edges(seq_edges, sv_edges)
    except Exception as e:
        logger.warning("Chromoauxesis classifier failed: {}".format(e))
        return {"decision": "not_chromoauxesis", "probability": 0.0, "features": {}}


def _get_intervals_from_seq_edges(seq_edges):
    intervals = defaultdict(list)
    for e in seq_edges:
        if e["chrom"]:
            intervals[e["chrom"]].append((e["start"], e["end"]))
    return dict(intervals)


def find_chromoauxesis_ecdna_cycles(cycleList, cycleCNs, segSeqD, graph_cns, invalidInds, bfb_cycle_inds):
    """Find circular cycles >250 kbp where every segment's graph CN exceeds the chromoauxesis ecDNA threshold."""
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
            if not seg_cns or max(seg_cns) <= ConfigVars.chromoauxesis_ecDNA_min_cn:
                all_above = False
                break
        if all_above:
            qualifying.append(ind)
    return qualifying


def get_raw_cycle_props(cycleList, maxCN, rearr_e, tot_over_min_cn, graph_cns):
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
        if cycleIsNoAmpInvalid(cycle, cycleCNs[ind], segSeqD, isSingleton, is_only_decomp, maxCN, graph_cns) and not args.force:
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


def run_classification(segSeqD, cycleList, cycleCNs, lc_filtered_cycle_count=0):
    graph_cns = get_graph_cns(graphFile, args.add_chr_tag)
    # first compute some properties about the foldbacks and copy numbers
    fb_edges, fb_readcount, fb_read_prop, maxCN, tot_over_min_cn = compute_f_from_AA_graph(graphFile)
    rearr_e = tot_rearr_edges(graphFile)
    (totalCompCyclicCont, totCyclicCont, ampClass, totalWeight, AMP_dvaluesDict, invalidInds, cycleTypes, cycleWeights,
     rearrCycleInds) = get_raw_cycle_props(cycleList, maxCN, rearr_e, tot_over_min_cn, graph_cns)

    # decomposition/amplicon complexity
    totalEnt, decompEnt, nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD, range(len(cycleList)),
                                                        set())
    AMP_dvaluesDict["Amp_entropy"] = totalEnt
    AMP_dvaluesDict["Amp_decomp_entropy"] = decompEnt
    AMP_dvaluesDict["Amp_nseg_entropy"] = nEnt

    fb_bwp, nfb_bwp, bfb_cwp, bfbHasEC, non_bfb_cycle_inds, bfb_cycle_inds = cycles_file_bfb_props(cycleList, segSeqD,
        cycleCNs, invalidInds, graphFile, graph_cns)

    # "foldback_read_prop", "BFB_bwp", "Distal_bwp", "BFB_cwp"
    AMP_dvaluesDict["foldback_read_prop"] = fb_read_prop
    AMP_dvaluesDict["BFB_bwp"] = fb_bwp
    AMP_dvaluesDict["Distal_bwp"] = nfb_bwp
    AMP_dvaluesDict["BFB_cwp"] = bfb_cwp
    bfbClass = classifyBFB(fb_read_prop, fb_bwp, nfb_bwp, bfb_cwp, maxCN, tot_over_min_cn)

    non_fb_rearr_e = rearr_e - fb_edges
    # heuristics to catch sequencing artifact samples
    if fb_edges > 15 and fb_read_prop > 0.8:
        bfbClass = False
        if non_fb_rearr_e >= 4 and tot_over_min_cn > ConfigVars.compCycContCut and maxCN > ConfigVars.sig_amp:
            ampClass = "Complex-non-cyclic"
        elif tot_over_min_cn > ConfigVars.compCycContCut and maxCN > ConfigVars.sig_amp:
            ampClass = "Linear"
        else:
            ampClass = INVALID_CLASS
            logger.warning(
                "{} has high foldback edge count ({}) and foldback read fraction ({:.3f}), "
                "suggesting read-orientation artifacts. Marking amplicon Invalid. Re-run the sample with "
                "AmpliconSuite-pipeline --foldback_pair_support_min 3 or higher if artifacts persist.".format(
                    amplicon_log_id(sName, ampN), fb_edges, fb_read_prop
                )
            )

    pre_finalize_ampClass = ampClass
    ampClass = finalize_no_amp_invalid_class(ampClass, lc_filtered_cycle_count)
    if pre_finalize_ampClass == NO_AMP_INVALID_INTERNAL_CLASS and ampClass == INVALID_CLASS:
        logger.info(
            "{} marked Invalid because {} cycle(s) were removed by low-complexity filtering.".format(
                amplicon_log_id(sName, ampN), lc_filtered_cycle_count
            )
        )

    bfbarchitect_summaries.append(run_bfbarchitect(graphFile, fb_edges, "{}_amplicon{}".format(sName, ampN)))

    # Chromoauxesis classification (graph already read once by parse_graph_edges_raw)
    _ca_seq_edges, _ca_sv_edges = parse_graph_edges_raw(graphFile, args.add_chr_tag)
    ca_result = _run_chromoauxesis_from_edges(_ca_seq_edges, _ca_sv_edges)
    ca_result["amplicon_intervals"] = _get_intervals_from_seq_edges(_ca_seq_edges)
    chromoauxesis_results.append(ca_result)
    is_chromoauxesis = (ca_result["decision"] == "chromoauxesis")

    ecStat = False
    bfbStat = False
    qualifying_ecdna_inds = []

    if is_chromoauxesis:
        # BFB detection is unchanged; ecDNA requires strict high-CN criterion
        if bfbClass and not is_non_feature_amplicon_class(ampClass):
            bfbStat = True
        else:
            bfb_cycle_inds = []
        if not is_non_feature_amplicon_class(ampClass):
            qualifying_ecdna_inds = find_chromoauxesis_ecdna_cycles(
                cycleList, cycleCNs, segSeqD, graph_cns, invalidInds, bfb_cycle_inds
            )
            if qualifying_ecdna_inds:
                ecStat = True

    else:
        if ampClass == "Cyclic" and not bfbClass:
            ecStat = True
            bfb_cycle_inds = []

        elif bfbClass and not is_non_feature_amplicon_class(ampClass):
            bfbStat = True
            if bfbHasEC:
                ecStat = True

        else:
            bfb_cycle_inds = []

    # determine number of ecDNA present (excluding BFB cycles)
    ecIndexClusters = []
    if ecStat:
        if is_chromoauxesis:
            non_qualifying = set(range(len(cycleList))) - set(qualifying_ecdna_inds)
            excludableCycleIndices = non_qualifying | set(bfb_cycle_inds) | set(invalidInds)
        else:
            excludableCycleIndices = set(bfb_cycle_inds + invalidInds)
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

    samp_to_ec_count[sName] += ecAmpliconCount
    # write entropy for each feature
    ecEntropies = []
    if ecAmpliconCount == 1 and not ecIndexClusters:
        ecEntropies.append((totalEnt, decompEnt, nEnt))

    for ecCycleList in ecIndexClusters:
        c_ex_I = bfb_cycle_inds if bfbStat else set()
        totalEnt, decompEnt, nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD, ecCycleList, c_ex_I)
        ecEntropies.append((totalEnt, decompEnt, nEnt))

    for ind, etup in enumerate(ecEntropies):
        featEntropyD[(sName, ampN, "ecDNA_" + str(ind + 1))] = etup

    if bfbStat:
        bfb_totalEnt, bfb_decompEnt, bfb_nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
                                                                        bfb_cycle_inds, set())
        featEntropyD[(sName, ampN, "BFB_1")] = (bfb_totalEnt, bfb_decompEnt, bfb_nEnt)

    if is_chromoauxesis:
        chromoauxesis_cycle_inds = set(range(len(cycleList))) - set(invalidInds)
        ca_totalEnt, ca_decompEnt, ca_nEnt = decompositionComplexity(
            graphFile, cycleList, cycleCNs, segSeqD, chromoauxesis_cycle_inds, set()
        )
        featEntropyD[(sName, ampN, "chromoauxesis_1")] = (ca_totalEnt, ca_decompEnt, ca_nEnt)

    bpg_linelist, gseg_cn_d, other_class_c_inds, feature_dict, prop_dict, ampClass = amplicon_annotation(cycleList, segSeqD,
        bfb_cycle_inds, ecIndexClusters, invalidInds, bfbStat, ecStat, ampClass, graphFile, args.add_chr_tag, lcD,
        args.ref, ca_result.get("amplicon_intervals") if is_chromoauxesis else None)

    bpgi_list.append(bpg_linelist)
    fd_list.append(feature_dict)
    prop_list.append(prop_dict)
    trim_sname = sName.rsplit("/")[-1].rsplit("_amplicon")[0]
    for feat_name, curr_fd in feature_dict.items():
        if curr_fd:
            full_fname = trim_sname + "_" + ampN + "_" + feat_name
            full_featname_to_graph[full_fname] = graphFile
            full_featname_to_intervals[full_fname] = curr_fd
            register_feature(feature_registry, full_fname, trim_sname, ampN, feat_name, graphFile, cyclesFile, curr_fd)

    if not bfbStat and not ecStat and not is_non_feature_amplicon_class(ampClass):
        featEntropyD[(sName, ampN, ampClass + "_1")] = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
            other_class_c_inds, set())

    feat_to_amped_genes = get_genes_from_intervals(gene_lookup, feature_dict, gseg_cn_d)
    feat_to_amped_ncrna = get_genes_from_intervals(ncRNA_lookup, feature_dict, gseg_cn_d)
    ftgd_list.append([sName, ampN, feat_to_amped_genes, feat_to_amped_ncrna])

    # store this additional information
    AMP_classifications.append((ampClass, ecStat, bfbStat, ecAmpliconCount))
    dvalues = [AMP_dvaluesDict[x] for x in categories]
    AMP_dvaluesList.append(dvalues)

    annotated_cycle_outname = os.path.basename(cyclesFile).rsplit("_cycles")[0] + "_annotated_cycles.txt"
    ftci_list.append([annotated_cycle_outname, cycleList, cycleCNs, segSeqD, bfb_cycle_inds, ecIndexClusters,
                      invalidInds, rearrCycleInds])


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
              "foldback_read_prop", "BFB_bwp", "Distal_bwp", "BFB_cwp", "Amp_entropy", "Amp_decomp_entropy",
              "Amp_nseg_entropy"]

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
    parser.add_argument("--plotstyle", help="Type of visualizations to produce.",
                        choices=["grouped", "individual", "noplot"], default="noplot")
    parser.add_argument("--force", help="Disable No-FSCNA/Invalid class if possible", action='store_true')
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files.",
                        action='store_true')
    parser.add_argument("--report_complexity", help="[Deprecated - on by default] Compute a measure of amplicon entropy"
                        " for each amplicon.", action='store_true', default=True)
    parser.add_argument("--verbose_classification", help="Generate verbose output with raw classification scores.",
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
    parser.add_argument("--make_results_table", help="Create summary results table after classification completes. "
                        "Only works when using --AA_results or --input (not with -c/-g single amplicon mode).", action='store_true')
    parser.add_argument("--bfbarchitect", help="Enable BFBArchitect score annotation. Results are reported in verbose "
                        "classification profiles and do not change classifications.", action='store_true')
    parser.add_argument("--bfb_threads", help="Number of threads for BFBArchitect ILP solver.", type=int, default=1)
    parser.add_argument("-v", "--version", action='version', version=__ampliconclassifier_version__)
    args = parser.parse_args()

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
    logger = setup_logger(args.o)

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
        logger.error(
            "--make_results_table can only be used with --AA_results or --input, not with single amplicon mode (-c/-g)")
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
        cmd = "{}/make_input.sh {} {}".format(src_dir, args.AA_results, outdir_loc + "/AC")
        logger.info("Generating .input file...")
        logger.info(cmd)

        # Capture stdout and stderr
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

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

    if args.bfbarchitect and args.ref == "hg19" and args.add_chr_tag:
        logger.error("--bfbarchitect cannot be used with --ref hg19 --add_chr_tag because BFBArchitect reads raw "
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

    if args.bfbarchitect:
        try:
            from bfbarchitect import reconstruct_bfb_from_graph, write_bfb_graph, write_bfb_cycles, visualize_BFB
            BFBARCHITECT_RECONSTRUCT = reconstruct_bfb_from_graph
            BFBARCHITECT_WRITE_GRAPH = write_bfb_graph
            BFBARCHITECT_WRITE_CYCLES = write_bfb_cycles
            BFBARCHITECT_VISUALIZE = visualize_BFB
            BFBARCHITECT_CENTROMERES = load_bfbarchitect_centromeres(AA_DATA_REPO_BASE, args.ref)
            if BFBARCHITECT_CENTROMERES:
                logger.info("BFBArchitect integration enabled.")

        except ImportError:
            logger.warning("BFBArchitect is not installed. Skipping BFBArchitect annotation.")

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

    # GLOBAL STORAGE VARIABLES
    ftgd_list = []  # store list of feature gene classifications
    ftci_list = []  # store list of cycles file info
    bpgi_list = []  # store list of bpg
    fd_list = []  # store list of feature_dicts
    prop_list = [] # store list of basic amplicon properties
    AMP_dvaluesList = []
    # EDGE_dvaluesList = []
    AMP_classifications = []
    bfbarchitect_summaries = []
    chromoauxesis_results = []
    sampNames = []
    cyclesFiles = []
    featEntropyD = {}
    samp_to_ec_count = defaultdict(int)
    samp_amp_to_graph = {}
    full_featname_to_graph = {}
    full_featname_to_intervals = {}
    feature_registry = {}

    # ITERATE OVER FILES CONDUCT THE CLASSIFICATION
    for i, fpair in enumerate(flist, 1):
        if len(fpair) > 2:
            orig_sName, cyclesFile, graphFile = fpair[:3]
            #sName = orig_sName.rsplit("_amplicon")[0]
            sName = re.split("_amplicon[0-9]*", orig_sName)[0]
            sampNames.append(sName)
            cyclesFiles.append(cyclesFile)
            ampN = cyclesFile.rstrip("_cycles.txt").rsplit("_")[-1]
            samp_amp_to_graph[sName + "_" + ampN] = graphFile
            logger.info("[{}/{}] {} {}".format(i, len(flist), sName, ampN))
            segSeqD, cycleList, cycleCNs, lc_filtered_cycle_count = parseCycle(cyclesFile, graphFile, add_chr_tag, lcD,
                                                                               patch_links)
            if args.ref == "hg19" or args.ref == "GRCh37":
                chroms = set([x[0] for k,x in segSeqD.items() if k != 0])
                if args.ref == "hg19" and not any([x.startswith("chr") for x in chroms]):
                    logger.warning("chrom names: " + str(chroms))
                    logger.warning("Warning! --ref hg19 was set, but chromosome naming is not consistent with hg19 ('chr1', 'chr2', etc.)\n")
                elif args.ref == "GRCh37" and any([x.startswith("chr") for x in chroms]) and not add_chr_tag:
                    logger.warning("chrom names: " + str(chroms))
                    logger.warning("Warning! --ref GRCh37 was set, but chromosome naming is not consistent with GRCh37 ('1', '2', etc.)\n")

        else:
            logger.info(str(fpair))
            logger.error("File list not properly formatted\n")
            sys.exit(1)

        run_classification(segSeqD, cycleList, cycleCNs, lc_filtered_cycle_count)

    logger.info("Classification stage completed")
    logger.info("Computing feature similarity scores...")
    feature_similarity_rows = compute_feature_similarity_scores(feature_registry)
    write_feature_similarity_scores(feature_similarity_rows, args.o)
    if args.filter_similar:
        logger.info("Filtering similar amplicons...")
        filter_similar_amplicons(len(flist), args.filter_pval, feature_registry, feature_similarity_rows)

    # make any requested visualizations
    plotting()

    logger.info("Writing outputs...")
    # OUTPUT FILE WRITING
    write_outputs(args, ftgd_list, ftci_list, bpgi_list, featEntropyD, categories, sampNames, cyclesFiles,
                  AMP_classifications, AMP_dvaluesList, samp_to_ec_count, fd_list, samp_amp_to_graph, prop_list,
                  summary_map, bfbarchitect_summaries, chromoauxesis_results)

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

    logger.info("")
    logger.info("Complete")
