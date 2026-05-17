#!/usr/bin/env python3

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"

import argparse
from collections import Counter
import contextlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import io
import logging
from math import log
import operator
import os
import re
import shlex
import subprocess
import sys

from intervaltree import IntervalTree, Interval

from ampclasslib.ac_annotation import *
from ampclasslib.ac_io import *
from ampclasslib.ac_util import ConfigVars, is_human_viral_hybrid
from ampclasslib.classification_records import (
    AmpliconInput, AmpliconRecord, ClassificationContext, Feature, collect_amplicon_records
)
from ampclasslib.config_params import load_config
from ampclasslib.radar_plotting import *
from ampclasslib._version import __ampliconclassifier_version__
from ampclasslib.chromoauxesis_ml import classify_from_edges as _chromoauxesis_classify_from_edges


add_chr_tag = False
WORKER_CONTEXT = None
BFBARCHITECT_RECONSTRUCT = None
BFBARCHITECT_WRITE_GRAPH = None
BFBARCHITECT_WRITE_CYCLES = None
BFBARCHITECT_VISUALIZE = None
BFBARCHITECT_CENTROMERES = None
NO_AMP_INVALID_INTERNAL_CLASS = "No amp/Invalid"
NO_FSCNA_CLASS = "No-FSCNA"
INVALID_CLASS = "Invalid"
BFBARCHITECT_SOURCE = "BFBArchitect"
AC_SOURCE = "AC"


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
def compute_f_from_AA_graph(graphf, lcD):
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
        "whole_graph_used": "NA",
        "passing_intervals": [],
        "passing_scores": [],
        "bfb_sources": "NA"
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
    passing_intervals = []
    passing_scores = []

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

    return {
        "min_score": format_bfbarchitect_score(min(all_scores) if all_scores else None),
        "passing_region_count": str(passing_regions),
        "multiplicities": ";".join(multiplicities) if multiplicities else "NA",
        "regions": ";".join(region_strings) if region_strings else "NA",
        "whole_graph_used": str(whole_graph_used),
        "passing_intervals": passing_intervals,
        "passing_scores": passing_scores,
        "bfb_sources": "NA"
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
            with suppress_bfbarchitect_chatter():
                BFBARCHITECT_WRITE_GRAPH(graph_out, res["new_segments"], res["svs"], res["sv_info"])
                BFBARCHITECT_WRITE_CYCLES(cycles_out, res["new_segments"], res["bfb_strings"], res["scores"],
                                          res["multiplicity"])
                BFBARCHITECT_VISUALIZE(cycle_file=cycles_out, graph_file=graph_out, cnr_file=None,
                                       output_prefix="{}_BFB".format(region_prefix), multiple=False,
                                       centromere_dict=BFBARCHITECT_CENTROMERES)
        except Exception as e:
            logger.warning("Could not write BFBArchitect outputs for {}: {}".format(region_prefix, str(e)))


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
        with suppress_bfbarchitect_chatter():
            results = BFBARCHITECT_RECONSTRUCT(graphf, centromere_dict=BFBARCHITECT_CENTROMERES, solver=None,
                                               multiple=False, whole_graph=False, threads=args.bfb_threads)
        whole_graph_used = False
        if should_retry_bfbarchitect_whole_graph(results):
            with suppress_bfbarchitect_chatter():
                results = BFBARCHITECT_RECONSTRUCT(graphf, centromere_dict=BFBARCHITECT_CENTROMERES, solver=None,
                                                   multiple=False, whole_graph=True, threads=args.bfb_threads)
            whole_graph_used = True

        summary = summarize_bfbarchitect_results(results, whole_graph_used)
        write_bfbarchitect_outputs(results, output_prefix, whole_graph_used)
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


def build_bfb_features(bfb_cycle_inds, bfbarchitect_summary, cycleList, segSeqD, invalidInds, graphFile, graph_cns):
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
        seq_edges, _sv_edges = parse_graph_edges_raw(graphFile, args.add_chr_tag)
        snapped_intervals = snap_bfbarchitect_intervals(passing_intervals, seq_edges, args.add_chr_tag)
        for ind, region in enumerate(snapped_intervals, 1):
            interval_dict = interval_dict_from_intervals([region])
            cycle_indices = cycle_indices_overlapping_bfb_intervals(
                cycleList, segSeqD, interval_dict, ConfigVars.bfbarchitect_cycle_overlap_threshold
            )
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


def plotting(classification_results):
    textCategories = ["No-FSCNA/\nInvalid", "Linear\namplification", "Trivial\ncycle", "Complex\nnon-cyclic",
                      "Complex\ncyclic", "BFB\nfoldback"]
    if args.plotstyle == "grouped":
        logger.info("plotting")
        make_classification_radar(
            textCategories, classification_results.AMP_dvaluesList, args.o + "_amp_class",
            classification_results.sampNames
        )

    elif args.plotstyle == "individual":
        logger.info("plotting")
        for a, s in zip(classification_results.AMP_dvaluesList, classification_results.sampNames):
            # print(textCategories, a)
            make_classification_radar(textCategories, [a[:len(textCategories)], ], args.o + "_" + s + "_amp_class",
                                      classification_results.sampNames)


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
    The following ClassificationResults fields are updated during similarity filtering:
    ftgd_list, ftci_list, bpgi_list, fd_list, prop_list, featEntropyD,
    AMP_classifications, chromoauxesis_results, and samp_to_ec_count.
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
            del classification_results.featEntropyD[(sname, ampN, fname + "_" + fnum)]

        except KeyError:
            pass

    for sname, ampN in amps_to_filter:
        keys_to_del = [x for x in classification_results.featEntropyD if x[0] == sname and x[1] == ampN]
        for key in keys_to_del:
            del classification_results.featEntropyD[key]

    for ind, sname in enumerate(classification_results.sampNames):
        ampN = classification_results.cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
        samp_amp = sname + "_" + ampN
        if samp_amp in samp_amp_filt_set:
            ampClass, ecStat, bfbStat, ecAmpliconCount = classification_results.AMP_classifications[ind]
            if (sname, ampN) in amps_to_filter:
                if ecAmpliconCount:
                    classification_results.samp_to_ec_count[sname] -= ecAmpliconCount

                classification_results.AMP_classifications[ind] = (NO_FSCNA_CLASS, False, False, 0)
                if ind < len(classification_results.chromoauxesis_results):
                    classification_results.chromoauxesis_results[ind]["decision"] = "not_chromoauxesis"
                    classification_results.chromoauxesis_results[ind]["probability"] = 0.0
                    classification_results.chromoauxesis_results[ind]["filtered_by_similarity"] = True
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


def get_raw_cycle_props(cycleList, cycleCNs, segSeqD, maxCN, rearr_e, tot_over_min_cn, graph_cns):
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
    global WORKER_CONTEXT, args, add_chr_tag, decomposition_strictness
    WORKER_CONTEXT = context
    args = context.args
    add_chr_tag = context.add_chr_tag
    decomposition_strictness = context.decomposition_strictness


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

    feature_entropy = {}
    feature_registry = {}
    full_featname_to_graph = {}
    full_featname_to_intervals = {}

    graph_cns = get_graph_cns(graphFile, args.add_chr_tag)
    # first compute some properties about the foldbacks and copy numbers
    fb_edges, fb_readcount, fb_read_prop, maxCN, tot_over_min_cn = compute_f_from_AA_graph(graphFile, lcD)
    rearr_e = tot_rearr_edges(graphFile)
    (totalCompCyclicCont, totCyclicCont, ampClass, totalWeight, AMP_dvaluesDict, invalidInds, cycleTypes, cycleWeights,
     rearrCycleInds) = get_raw_cycle_props(cycleList, cycleCNs, segSeqD, maxCN, rearr_e, tot_over_min_cn, graph_cns)

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
    bfb_suppressed_as_artifact = False
    # heuristics to catch sequencing artifact samples
    if fb_edges > 15 and fb_read_prop > 0.8:
        bfbClass = False
        bfb_suppressed_as_artifact = True
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

    if bfb_suppressed_as_artifact:
        bfbarchitect_summary = empty_bfbarchitect_summary()
    else:
        bfbarchitect_summary = run_bfbarchitect(graphFile, fb_edges, "{}_{}".format(sName, ampN))

    native_bfb_cycle_inds = list(bfb_cycle_inds) if bfbClass else []
    bfb_features = build_bfb_features(
        native_bfb_cycle_inds, bfbarchitect_summary, cycleList, segSeqD, invalidInds, graphFile, graph_cns
    )
    bfb_cycle_inds = sorted(set().union(*[feature.cycle_indices for feature in bfb_features])) if bfb_features else []
    bfb_sources = sorted(set().union(*[feature.sources for feature in bfb_features])) if bfb_features else []
    if bfbClass and AC_SOURCE not in bfb_sources:
        bfb_sources.append(AC_SOURCE)
    bfbarchitect_summary["bfb_sources"] = "|".join(sorted(bfb_sources)) if bfb_sources else "NA"
    bfb_detected = bool(bfbClass) or bool(bfb_features)
    bfbarchitect_detected = BFBARCHITECT_SOURCE in bfb_sources

    # Chromoauxesis classification (graph already read once by parse_graph_edges_raw)
    _ca_seq_edges, _ca_sv_edges = parse_graph_edges_raw(graphFile, args.add_chr_tag)
    ca_result = _run_chromoauxesis_from_edges(_ca_seq_edges, _ca_sv_edges)
    ca_result["amplicon_intervals"] = _get_intervals_from_seq_edges(_ca_seq_edges)
    is_chromoauxesis = (ca_result["decision"] == "chromoauxesis")

    ecStat = False
    bfbStat = False
    qualifying_ecdna_inds = []

    if is_chromoauxesis:
        # BFB detection is unchanged; ecDNA requires strict high-CN criterion
        if bfb_detected and not is_non_feature_amplicon_class(ampClass):
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
        if ampClass == "Cyclic" and not bfb_detected:
            ecStat = True
            bfb_cycle_inds = []

        elif bfb_detected and not is_non_feature_amplicon_class(ampClass):
            bfbStat = True
            if bfbHasEC or (bfbarchitect_detected and ampClass == "Cyclic"):
                ecStat = True

        else:
            bfb_cycle_inds = []
            bfb_features = []

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
        feature_entropy[(sName, ampN, "ecDNA_" + str(ind + 1))] = etup

    if bfbStat:
        bfb_entropy_cycle_inds = bfb_cycle_inds if bfb_cycle_inds else range(len(cycleList))
        bfb_totalEnt, bfb_decompEnt, bfb_nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
                                                                        bfb_entropy_cycle_inds, set())
        feature_entropy[(sName, ampN, "BFB_1")] = (bfb_totalEnt, bfb_decompEnt, bfb_nEnt)

    if is_chromoauxesis:
        chromoauxesis_cycle_inds = set(range(len(cycleList))) - set(invalidInds)
        ca_totalEnt, ca_decompEnt, ca_nEnt = decompositionComplexity(
            graphFile, cycleList, cycleCNs, segSeqD, chromoauxesis_cycle_inds, set()
        )
        feature_entropy[(sName, ampN, "chromoauxesis_1")] = (ca_totalEnt, ca_decompEnt, ca_nEnt)

    bpg_linelist, gseg_cn_d, other_class_c_inds, feature_dict, prop_dict, ampClass = amplicon_annotation(cycleList, segSeqD,
        bfb_cycle_inds, ecIndexClusters, invalidInds, bfbStat, ecStat, ampClass, graphFile, args.add_chr_tag, lcD,
        args.ref, ca_result.get("amplicon_intervals") if is_chromoauxesis else None, bfb_features)

    trim_sname = sName.rsplit("/")[-1].rsplit("_amplicon")[0]
    for feat_name, curr_fd in feature_dict.items():
        if curr_fd:
            full_fname = trim_sname + "_" + ampN + "_" + feat_name
            full_featname_to_graph[full_fname] = graphFile
            full_featname_to_intervals[full_fname] = curr_fd
            register_feature(feature_registry, full_fname, trim_sname, ampN, feat_name, graphFile, cyclesFile, curr_fd)

    if not bfbStat and not ecStat and not is_non_feature_amplicon_class(ampClass):
        feature_entropy[(sName, ampN, ampClass + "_1")] = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
            other_class_c_inds, set())

    feat_to_amped_genes = get_genes_from_intervals(gene_lookup, feature_dict, gseg_cn_d)
    feat_to_amped_ncrna = get_genes_from_intervals(ncRNA_lookup, feature_dict, gseg_cn_d)
    feature_gene_entry = [sName, ampN, feat_to_amped_genes, feat_to_amped_ncrna]

    # store this additional information
    classification = (ampClass, ecStat, bfbStat, ecAmpliconCount)
    dvalues = [AMP_dvaluesDict[x] for x in categories]

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
        chromoauxesis_result=ca_result,
        feature_entropy=feature_entropy,
        ec_count=ecAmpliconCount,
        samp_amp_to_graph={sName + "_" + ampN: graphFile},
        feature_registry=feature_registry,
        full_featname_to_graph=full_featname_to_graph,
        full_featname_to_intervals=full_featname_to_intervals,
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
    parser.add_argument("--no_bfbarchitect", help="Disable BFBArchitect integration. By default, AmpliconClassifier "
                        "uses BFBArchitect when it is installed and available.", action='store_true')
    parser.add_argument("--bfb_threads", help="Number of threads for BFBArchitect ILP solver.", type=int, default=1)
    parser.add_argument("--jobs", help="Number of amplicons to classify in parallel.", type=int, default=1)
    parser.add_argument("-v", "--version", action='version', version=__ampliconclassifier_version__)
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
            from bfbarchitect import reconstruct_bfb_from_graph, write_bfb_graph, write_bfb_cycles, visualize_BFB
            BFBARCHITECT_RECONSTRUCT = reconstruct_bfb_from_graph
            BFBARCHITECT_WRITE_GRAPH = write_bfb_graph
            BFBARCHITECT_WRITE_CYCLES = write_bfb_cycles
            BFBARCHITECT_VISUALIZE = visualize_BFB
            BFBARCHITECT_CENTROMERES = load_bfbarchitect_centromeres(AA_DATA_REPO_BASE, args.ref)
            if BFBARCHITECT_CENTROMERES:
                args.bfbarchitect = True
                logger.info("BFBArchitect integration enabled.")
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
    classification_results = collect_amplicon_records(amplicon_records)

    logger.info("Computing feature similarity scores...")
    feature_similarity_rows = compute_feature_similarity_scores(classification_results.feature_registry)
    write_feature_similarity_scores(feature_similarity_rows, args.o)
    if args.filter_similar:
        logger.info("Filtering similar amplicons...")
        filter_similar_amplicons(
            len(flist), args.filter_pval, classification_results.feature_registry,
            feature_similarity_rows, classification_results
        )

    # make any requested visualizations
    plotting(classification_results)

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

    logger.info("")
    logger.info("Complete")
