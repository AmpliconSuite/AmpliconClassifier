#!/usr/bin/env python3

__version__ = "0.3.5"
__author__ = "Jens Luebeck"

import argparse
from collections import defaultdict
import copy
from math import log
import operator
import os
import sys

from intervaltree import IntervalTree

tot_min_del = 5000  # minimum size of deletion before non-trivial
minCycleSize = 10000
compCycContCut = 50000
cycCut = 0.12
compCut = 0.3
min_upper_cn = 4.5

# bfb thresholds
min_score_for_bfb = 0.25
min_fb_reads_for_bfb = 10
fb_dist_cut = 25000


class breakpoint(object):
    def __init__(self, lchrom, lpos, rchrom, rpos, cn):
        self.lchrom = lchrom
        self.lpos = lpos
        self.rchrom = rchrom
        self.rpos = rpos
        self.cn = cn

    def to_string(self):
        return self.lchrom + ":" + str(self.lpos) + " | " + self.rchrom + ":" + str(self.rpos) + "\t" + str(self.cn)


# ------------------------------------------------------------
# Methods to compute values used in classification
def get_size(cycle, segSeqD):
    return sum(segSeqD[abs(x)][2] - segSeqD[abs(x)][1] for x in cycle)


def weightedCycleAmount(cycle, cn, segSeqD):
    # get length of cycle
    sc_length = get_size(cycle, segSeqD) / 1000.
    return sc_length * cn


def get_diff(e1, e2, segSeqD):
    p1_abs = segSeqD[abs(e1)]
    p2_abs = segSeqD[abs(e2)]
    if e1 == 0 or e2 == 0:
        return 1

    p1_end = p1_abs[2] if e1 > 0 else p1_abs[1]
    p2_start = p2_abs[1] if e2 > 0 else p2_abs[2]
    return abs(p2_start - p1_end)


def isCircular(cycle):
    circular = False if cycle[0] == 0 and cycle[-1] == 0 else True
    return circular


def isRearranged(cycle, segSeqD):
    # check if it contains regions from multiple chroms
    chromList = [segSeqD[abs(ind)][0] for ind in cycle if ind != 0]
    if len(set(chromList)) > 1:
        return True

    max_del_size = 0
    for i in range(0, len(cycle) - 1):
        if cycle[i] == 0 or cycle[i + 1] == 0:
            continue
        if cycle[i] < 0 and cycle[i + 1] > 0 or cycle[i] > 0 and cycle[i + 1] < 0:
            return True

        dist_diff = get_diff(cycle[i], cycle[i + 1], segSeqD)
        max_del_size = max(dist_diff, max_del_size)
        # tot_del_size += dist_diff
        # if tot_del_size > tot_min_del:
        #     print("Delsize")
        #     return True
    if max_del_size > tot_min_del:
        return True

    return False


def tot_rearr_edges(graphf, add_chr_tag):
    rearr_e = 0
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                lbp, rbp = fields[1].split("->")
                lchrom, lpd = lbp.rsplit(":")
                rchrom, rpd = rbp.rsplit(":")
                if add_chr_tag and not lchrom.startswith("chr"):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                lpos, ldir = int(lpd[:-1]), lpd[-1]
                rpos, rdir = int(rpd[:-1]), rpd[-1]
                if ldir == rdir:
                    # if lchrom == rchrom and abs(rpos - lpos) < fb_dist_cut:
                    rearr_e += 1

                elif abs(rpos - lpos) > fb_dist_cut:
                    rearr_e += 1

    return rearr_e


def decompositionComplexity(graphf, cycleList, cycleCNs, segSeqD, feature_inds, exclude_inds):
    #construct intervaltree of valid regions
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
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("sequence"):
                fields = line.rsplit()
                c, s, e = fields[1].rsplit(":")[0], int(fields[1].rsplit(":")[1][:-1]), int(fields[2].rsplit(":")[1][:-1])+1
                if not hit_region_it[c][s:e]:
                    continue

                cn = float(fields[3])
                size = float(fields[5]) / 1000.

                # if cn > 1:
                segs += 1
                totalGraphWeight += (size * cn)

            elif line.startswith("BreakpointEdge"):
                break

    # cycleWeights = [0] * len(feature_inds)
    # for ind, cycle in enumerate(cycleList):
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
                wca = weightedCycleAmount(cycle, cycleCNs[ind], segSeqD)
                if ind in feature_inds:
                    new_feat_inds.add(len(cycleWeights))

                cycleWeights.append(wca)


    # scW = sorted(cycleWeights, reverse=True)
    # cf = cycleWeights[0]/totalGraphWeight
    cf = 0
    fe_ent = 0
    added_cf = 0
    cInd = 0
    while cf + added_cf < hf_cut and cInd < len(cycleWeights):
        if cInd in new_feat_inds:
            cf += added_cf
            if added_cf > 0:
                fe_ent += (added_cf * log(added_cf))

            added_cf = cycleWeights[cInd] / float(totalGraphWeight)

        cInd += 1

    cf+=added_cf
    cf = round(cf, 5)
    rf = (1 - cf)
    # print(rf, cf, totalGraphWeight)
    if rf > 0:
        fu_ent = -1 * rf * log(rf)
    else:
        fu_ent = 0

    seg_ent = log(1.0 / segs) if segs > 0 else 0
    return fu_ent - fe_ent - seg_ent, fu_ent - fe_ent, -1 * seg_ent


# Compute f (foldback fraction) from the edges in the AA graph alone
def compute_f_from_AA_graph(graphf, add_chr_tag):
    with open(graphf) as infile:
        fbCount, nonFbCount, fbEdges, maxCN = 0, 0, 0, 0
        for line in infile:
            fields = line.rstrip().rsplit()
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

                rSupp = int(fields[3])
                if ldir == rdir:
                    if lchrom == rchrom and abs(rpos - lpos) < fb_dist_cut:
                        fbCount += rSupp
                        fbEdges += 1

                    else:
                        nonFbCount += rSupp

                else:
                    nonFbCount += rSupp

            elif line.startswith("sequence"):
                if not lcD[fields[1].rsplit(":")[0]].overlaps(int(fields[1].rsplit(":")[1][:-1]),
                                                              int(fields[2].rsplit(":")[1][:-1])):

                    ccn = float(fields[3])
                    if ccn > maxCN:
                        maxCN = ccn

    # just return 0 if there isn't enough support
    if fbEdges < 2:
        return 0, maxCN

    return fbCount / max(1.0, fbCount + nonFbCount), maxCN


def nonbfb_cycles_are_ecdna(non_bfb_cycle_inds, cycleList, segSeqD, cycleCNs):
    for ind in non_bfb_cycle_inds:
        cycle = cycleList[ind]
        length = get_size(cycle, segSeqD)

        if length > 100000 and cycleCNs[ind] > 5:
            return True

    return False


# proportion of cycles with foldbacks
def cycles_file_bfb_props(cycleList, segSeqD, cycleCNs):
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

        hit = False
        isBFBelem = False
        illegalBFB = False
        for a, b in zip(cycle[:-1], cycle[1:]):
            # changes direction on same chrom
            diff = get_diff(a, b, segSeqD)
            aSize = get_size([a, ], segSeqD)
            bSize = get_size([b, ], segSeqD)
            if aSize < minCycleSize and bSize < minCycleSize:
                continue

            if a * b < 0 and segSeqD[abs(a)][0] == segSeqD[abs(b)][0]:
                hit = True
                if diff < 50000:
                    isBFBelem = True
                    FB_breaks += cycleCNs[ind]

                else:
                    distal_breaks += cycleCNs[ind]

            elif diff > tot_min_del:
                hit = True
                distal_breaks += cycleCNs[ind]

            if segSeqD[abs(a)][0] != segSeqD[abs(b)][0] and not (a == 0 or b == 0):
                illegalBFB = True

        if illegalBFB:
            isBFBelem = False

        if cycle[0] == 0 and not hit:
            lin_breaks += cycleCNs[ind]

        if isBFBelem:
            tot_bfb_supp_cycles += 1
            bfb_weight += cycleCNs[ind]
            bfb_cycle_inds.append(ind)

        elif cycle[0] != 0 and get_size(cycle[:-1], segSeqD) > 30000:
            non_bfb_cycle_weight += cycleCNs[ind]
            non_bfb_cycle_inds.append(ind)

    hasEC = nonbfb_cycles_are_ecdna(non_bfb_cycle_inds, cycleList, segSeqD, cycleCNs)
    # if len(cycleList) >= 5:
    #     minBFBCyclesRequired = 2
    # else:
    minBFBCyclesRequired = 2
    if FB_breaks > 1.5 and tot_bfb_supp_cycles >= minBFBCyclesRequired:
        tot = float(FB_breaks + distal_breaks + lin_breaks)
        return FB_breaks / tot, distal_breaks / tot, bfb_weight / (non_bfb_cycle_weight + bfb_weight), hasEC, \
               non_bfb_cycle_inds, bfb_cycle_inds

    return 0, 0, 0, False, [], []


# ------------------------------------------------------------
# Classifications
def cycleIsNoAmpInvalid(cycle, cn, segSeqD, isSingleton, maxCN):
    # CN flow can be split across multiple amps
    if not isSingleton:
        scale = min(args.min_cn_flow, maxCN / 10.)
    elif maxCN > 7:
        scale = min(3, maxCN / 8.)
    else:
        scale = 2.5

    if (cn <= scale) or (maxCN < min_upper_cn):
        return True

    length = get_size(cycle, segSeqD)
    return length < minCycleSize


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


# categories = ["No amp/Invalid", "Linear amplification", "Trivial cycle", "Complex non-cyclic", "Complex cyclic"]
def classifyAmpliconProfile(amp_profile, rearr_e, totalCompCycCont, force=False):
    cycSig = amp_profile["Trivial cycle"] + amp_profile["Complex cyclic"]
    if cycSig > cycCut or totalCompCyclicCont > compCycContCut:
        return "Cyclic"

    elif amp_profile["Complex non-cyclic"] + cycSig > compCut:
        if rearr_e < 2:
            return "Linear amplification"

        return "Complex non-cyclic"

    else:
        if max(amp_profile.values()) == 0:
            return "No amp/Invalid"

        elif amp_profile["No amp/Invalid"] > 0:
            if amp_profile["Linear amplification"] / float(amp_profile["No amp/Invalid"]) > 0.25:
                if rearr_e >= 5:
                    return "Complex non-cyclic"

                return "Linear amplification"

        if force:
            del amp_profile["No amp/Invalid"]
            if cycSig > max(amp_profile.values()):
                return "Cyclic"

        maxCat = max(amp_profile.items(), key=operator.itemgetter(1))[0]
        return maxCat


def classifyBFB(fb, cyc_sig, nonbfb_sig, bfb_cyc_ratio, maxCN):
    if fb < min_score_for_bfb or cyc_sig < 0.295 or maxCN < 4:
        return None

    # dominated by non-classical BFB cycles
    elif nonbfb_sig > 0.5 and bfb_cyc_ratio < 0.6:
        return None

    # if bfb_cyc_ratio < 0.85:
    #     # if nonbfb_sig > 0.55:
    #     #     return None
    #
    #     if args.use_BFB_linked_cyclic_class:
    #         return "BFB-linked cyclic"

    return "BFB"


# ------------------------------------------------------------
# # structure metanalysis
def clusterECCycles(cycleList, cycleCNs, segSeqD, excludableCycleIndices=None):
    padding = 500000
    indices = [x for x in range(len(cycleList)) if cycleList[x][0] != 0 and x not in excludableCycleIndices]
    clusters = []
    seenSegs = set()
    for ind in indices:
        cycle = cycleList[ind]
        if cycleCNs[ind] < args.min_cn_flow and get_size(cycle, segSeqD) < minCycleSize:
            continue

        cIndsToMerge = set()
        s_set = set([segSeqD[abs(s_num)] for s_num in cycle])
        s_set -= seenSegs
        if not s_set:
            continue

        for c_ind, clust_dict in enumerate(clusters):
            for s in s_set:
                if clust_dict[s[0]][s[1] - padding:s[2] + padding]:
                    cIndsToMerge.add(c_ind)
                    break

        # if cIndsToMerge:
        #     for s in s_set:
        #         clusters[cIndsToMerge
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
                currIndexSet.add(ival.data)  # set to 1-based

        indexClusters.append(currIndexSet)

    return indexClusters


# ------------------------------------------------------------
# Input and data prep
def bpgEdgeToCycles(bp, posCycleLookup):
    lCycles = set([x.data for x in posCycleLookup[bp.lchrom][bp.lpos]])
    rCycles = set([x.data for x in posCycleLookup[bp.rchrom][bp.rpos]])
    return lCycles, rCycles


# some have 0's in the middle because of an AA bug. This is fixes that automatically
def repair_cycle(cycle):
    if 0 in cycle and cycle[0] != 0:
        zero_ind = cycle.index(0)
        rotCyc = cycle[zero_ind + 1:] + cycle[:zero_ind + 1]

        return rotCyc

    return cycle


def parseCycle(cyclef, add_chr_tag):
    segSeqD = {0: (None, 0, 0)}
    cycleList = []
    cycleCNs = []
    seenCycs = set()

    with open(cyclef) as infile:
        for line in infile:
            if line.startswith("Segment"):
                fields = line.rstrip().rsplit()
                segNum = int(fields[1])
                chrom = fields[2]
                if add_chr_tag and not chrom.startswith('chr'):
                    chrom = "chr" + chrom

                l, r = int(fields[3]), int(fields[4])
                segSeqD[segNum] = (chrom, l, r)

            elif line.startswith("Cycle"):
                cf = [tuple(x.rsplit("=")) for x in line.rstrip().rsplit(";")]
                cd = dict(cf)
                ss = cd["Segments"]
                num_ss = [int(x[-1] + x[:-1]) for x in ss.rsplit(",")]
                # print(num_ss)
                # if any([segSeqD[abs(x)][0] == "hs37d5" for x in num_ss if x != 0]):
                #     continue
                lcCycle = False
                pop_inds = []
                for seg_ind, seg in enumerate(num_ss):
                    t = segSeqD[abs(seg)]
                    if lcD[t[0]].overlaps(t[1], t[2]):
                        if num_ss[0] == 0 and (seg_ind == 1 or seg_ind == len(num_ss) - 2) and len(num_ss) > 3:
                            pop_inds.append(seg_ind)
                            continue

                        else:
                            print("Cycle was LC", str(t[0]), str(t[1]), str(t[2]))
                            lcCycle = True
                            break

                if lcCycle:
                    continue

                elif pop_inds:
                    for seg_ind in pop_inds[::-1]:
                        num_ss.pop(seg_ind)


                currCycle = repair_cycle(num_ss)
                uid = ss + "," + cd["Copy_count"]
                if uid in seenCycs:
                    print(cyclef + " duplicate cycle encountered")

                else:
                    cycleList.append(currCycle)
                    seenCycs.add(uid)
                    cycleCNs.append(float(cd["Copy_count"]))

    return segSeqD, cycleList, cycleCNs


def parseBPG(bpgf, add_chr_tag):
    bps = []
    with open(bpgf) as infile:
        for line in infile:
            # if line.startswith("discordant") or line.startswith("concordant"):
            if line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                l, r = fields[1].rsplit("->")
                lchrom, lpos = l[:-1].rsplit(":")
                rchrom, rpos = r[:-1].rsplit(":")
                lpos, rpos = int(lpos), int(rpos)
                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                if lcD[lchrom][lpos] or lcD[rchrom][rpos]:
                    continue

                cn = float(fields[2])
                currBP = breakpoint(lchrom, lpos, rchrom, rpos, cn)
                bps.append(currBP)

    return bps


# build a lookup for position to list of cycles hitting it
def buildPosCycleLookup(cycles, segSeqD):
    posCycleLookup = defaultdict(IntervalTree)
    for ind, c in enumerate(cycles):
        for seg in c:
            if seg != 0:
                chrom, s, e = segSeqD[abs(seg)]
                posCycleLookup[chrom][s:e + 1] = ind

    return posCycleLookup


def buildLCDatabase(mappabilityFile):
    lcD = defaultdict(IntervalTree)
    with open(mappabilityFile) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            chrom, s, e = fields[0], int(fields[1]), int(fields[2])
            if e - s > 7500:
                lcD[chrom].addi(s, e)

    return lcD


def readFlist(filelist):
    flist = []
    with open(filelist) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                if len(fields) < 2 or len(fields) > 3:
                    print("Bad formatting in: ", line)
                else:
                    flist.append(fields)

    return flist


# ------------------------------------------------------------

'''
Amplicon Classes:
#if not invalid

1) No amp/Invalid
2) Linear amplification
3) Complex non-cyclic
4) BFB
5) BFB-linked cylic (ecBFB)
6) Cyclic (ecDNA)

Graph edge classes:
1) No amp/Invalid
2) Non-cyclic
3) Integration
4) Hybrid - joins amplified cyclic and non cyclic
5) Cyclic

'''

mixLookups = {
    frozenset(["No amp/Invalid", ]): "No amp/Invalid",
    frozenset(["No amp/Invalid", "Linear amplification"]): "Integration",
    frozenset(["No amp/Invalid", "Trivial cycle"]): "Integration",
    frozenset(["No amp/Invalid", "Complex non-cyclic"]): "Integration",
    frozenset(["No amp/Invalid", "Complex cyclic"]): "Integration",
    frozenset(["Linear amplification"]): "Non-cyclic",
    frozenset(["Linear amplification", "Trivial cycle"]): "Integration",
    frozenset(["Linear amplification", "Complex non-cyclic"]): "Non-cyclic",
    frozenset(["Linear amplification", "Complex cyclic"]): "Integration",
    frozenset(["Trivial cycle"]): "Cyclic",
    frozenset(["Trivial cycle", "Complex non-cyclic"]): "Hybrid",
    frozenset(["Trivial cycle", "Complex cyclic"]): "Cyclic",
    frozenset(["Complex non-cyclic"]): "Non-cyclic",
    frozenset(["Complex non-cyclic", "Complex cyclic"]): "Hybrid",
    frozenset(["Complex cyclic"]): "Cyclic",
}

# (circular,complex)
categories = ["No amp/Invalid", "Linear amplification", "Trivial cycle", "Complex non-cyclic", "Complex cyclic",
              "foldback_read_prop", "BFB_bwp", "Distal_bwp", "BFB_cwp", "Amp_entropy", "Amp_decomp_entropy",
              "Amp_nseg_entropy"]
mixing_cats = ["No amp/Invalid", "Non-cyclic", "Integration", "Hybrid", "Cyclic"]

ampDefs = {(False, False): "Linear amplification", (False, True): "Complex non-cyclic",
           (True, False): "Trivial cycle", (True, True): "Complex cyclic"}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify AA amplicon type")
    parser.add_argument("-c", "--cycles", help="AA-formatted cycles file")
    parser.add_argument("-g", "--graph", help="AA-formatted graph file")
    parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, or GRCh38",
                        choices=["hg19", "GRCh37", "hg38", "GRCh38"], required=True)

    parser.add_argument("--min_cn_flow", type=float, help="Minimum CN flow to consider as amplification", default=1)
    parser.add_argument("--min_size", type=float, help="Minimum cycle size (in bp) to consider as valid amplicon",
                        default=5000)
    parser.add_argument("-o", help="Output filename prefix")
    parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: \
    sample_name cycles.txt graph.txt")
    parser.add_argument("--plotStyle", help="Type of visualizations to produce",
                        choices=["grouped", "individual", "noplot"], default="noplot")
    parser.add_argument("--force", help="Disable No amp/Invalid class if possible", action='store_true')
    # parser.add_argument("--use_BFB_linked_cyclic_class", help="Include the \'BFB-linked cyclic\' class",
    #                     action='store_true')
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files",
                        action='store_true')
    parser.add_argument("--report_genes", help="Extract list of genes from amplicons with given classification.",
                        choices=["ecdna", "bfb", "other", "all"], nargs='+', default=[])
    parser.add_argument("--report_complexity", help="Compute a measure of amplicon entropy for each amplicon.",
                        action='store_true')
    parser.add_argument("--verbose_classification", help="Generate verbose output with raw classification scores.",
                        action='store_true')
    parser.add_argument("-v", "--version", action='version', version='amplicon_classifier {version} \n Author: Jens \
                        Luebeck (jluebeck [at] ucsd.edu)'.format(version=__version__))

    args = parser.parse_args()

    if not (args.cycles and args.graph) and not args.input:
        print("Need to specify (--cycles & --graph) or --input\n")
        sys.exit(1)

    if args.ref == "hg38": args.ref = "GRCh38"

    # check if aa data repo set, construct LC datatabase
    try:
        AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + args.ref + "/"
        fDict = {}
        with open(AA_DATA_REPO + "file_list.txt") as infile:
            for line in infile:
                fields = line.strip().rsplit()
                fDict[fields[0]] = fields[1]

        lcPath = AA_DATA_REPO + fDict["mapability_exclude_filename"]
        lcD = buildLCDatabase(lcPath)

    except KeyError:
        sys.stderr.write("$AA_DATA_REPO not set. Please see AA installation instructions.\n")
        sys.exit(1)

    gene_lookup = {}
    ftgd_list = []
    if args.report_genes:
        # read the gene list and import
        import get_genes
        gene_file_location_lookup = {"hg19": "human_hg19_september_2011/Genes_July_2010_hg19.gff",
                                     "GRCh38": "genes_hg38.gff",
                                     "GRCh37": "human_hg19_september_2011/Genes_July_2010_hg19.gff"}

        refGeneFileLoc = AA_DATA_REPO + gene_file_location_lookup[args.ref]
        gene_lookup = get_genes.parse_genes(refGeneFileLoc)

    if not args.input:
        tempName = args.cycles.rsplit("/")[-1].rsplit(".")[0]
        flist = [[tempName, args.cycles, args.graph]]
        if not args.o:
            args.o = os.path.basename(args.cycles).rsplit("_cycles.txt")[0]

    else:
        flist = readFlist(args.input)
        if not args.o:
            args.o = os.path.basename(args.input).rsplit(".")[0]


    minCycleSize = args.min_size

    AMP_dvaluesList = []
    EDGE_dvaluesList = []
    AMP_classifications = []
    sampNames = []
    cyclesFiles = []
    featEntropyD = {}
    for fpair in flist:
        if len(fpair) > 2:
            sName, cyclesFile, graphFile = fpair
            sampNames.append(sName)
            cyclesFiles.append(cyclesFile)
            segSeqD, cycleList, cycleCNs = parseCycle(cyclesFile, args.add_chr_tag)

        else:
            print(fpair)
            sys.stderr.write("File list not properly formatted\n")
            sys.exit(1)

        ampN = cyclesFile.rstrip("_cycles.txt").rsplit("_")[-1]
        print(sName, ampN)
        cycleTypes = []
        cycleWeights = []
        invalidInds = []
        fb_prop, maxCN = compute_f_from_AA_graph(graphFile, args.add_chr_tag)
        rearr_e = tot_rearr_edges(graphFile, args.add_chr_tag)
        totalCompCyclicCont = 0
        for ind, cycle in enumerate(cycleList):
            hasNonCircLen1 = True if len(cycle) == 3 and cycle[0] == 0 else False
            oneCycle = (len(cycleList) == 1)
            isSingleton = hasNonCircLen1 or oneCycle
            if cycleIsNoAmpInvalid(cycle, cycleCNs[ind], segSeqD, isSingleton, maxCN) and not args.force:
                invalidInds.append(ind)
                cycleTypes.append("No amp/Invalid")

            else:
                circCyc = isCircular(cycle)
                compCyc = isRearranged(cycle, segSeqD)
                if circCyc and compCyc:
                    totalCompCyclicCont += get_size(cycle, segSeqD)

                cycleTypes.append(ampDefs[(circCyc, compCyc)])

            currWt = weightedCycleAmount(cycle, cycleCNs[ind], segSeqD)
            cycleWeights.append(currWt)

        totalWeight = max(sum(cycleWeights), 1)
        AMP_dvaluesDict = {x: 0.0 for x in categories}
        for i, wt in zip(cycleTypes, cycleWeights):
            AMP_dvaluesDict[i] += (wt / totalWeight)

        # anything stored in AMP_dvaluesDict prior to running classify will get used in classification
        # make sure you're not putting in other properties before here.
        ampClass = classifyAmpliconProfile(AMP_dvaluesDict, rearr_e, totalCompCyclicCont)

        # decomposition/amplicon complexity
        totalEnt, decompEnt, nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
                                                                range(len(cycleList)), set())
        AMP_dvaluesDict["Amp_entropy"] = totalEnt
        AMP_dvaluesDict["Amp_decomp_entropy"] = decompEnt
        AMP_dvaluesDict["Amp_nseg_entropy"] = nEnt

        # now layer on the bfb classification
        # first compute some properties
        fb_prop, maxCN = compute_f_from_AA_graph(graphFile, args.add_chr_tag)

        fb_bwp, nfb_bwp, bfb_cwp, bfbHasEC, non_bfb_cycle_inds, \
        bfb_cycle_inds = cycles_file_bfb_props(cycleList, segSeqD, cycleCNs)
        # "foldback_read_prop", "BFB_bwp", "Distal_bwp", "BFB_cwp"
        AMP_dvaluesDict["foldback_read_prop"] = fb_prop
        AMP_dvaluesDict["BFB_bwp"] = fb_bwp
        AMP_dvaluesDict["Distal_bwp"] = nfb_bwp
        AMP_dvaluesDict["BFB_cwp"] = bfb_cwp

        bfbClass = classifyBFB(fb_prop, fb_bwp, nfb_bwp, bfb_cwp, maxCN)

        ecStat = False
        bfbStat = False
        if ampClass == "Cyclic" and not bfbClass:
            ecStat = True
            bfb_cycle_inds = []

        elif bfbClass:
            bfbStat = True
            if bfbHasEC:
                ecStat = True

        else:
            bfb_cycle_inds = []

        # determine number of ecDNA present
        ecIndexClusters = []
        if ecStat:
            excludableCycleIndices = set(bfb_cycle_inds + invalidInds)
            ecIndexClusters = clusterECCycles(cycleList, cycleCNs, segSeqD, excludableCycleIndices)
            ecAmpliconCount = max(len(ecIndexClusters), 1)

        else:
            ecAmpliconCount = 0

        # write entropy for each feature
        ecEntropies = []
        if ecAmpliconCount == 1 and not ecIndexClusters:
            ecEntropies.append((totalEnt, decompEnt, nEnt))

        for ecCycleList in ecIndexClusters:
            c_ex_I = bfb_cycle_inds if bfbStat else set()
            totalEnt, decompEnt, nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
                                                                ecCycleList, c_ex_I)
            ecEntropies.append((totalEnt, decompEnt, nEnt))

        for ind, etup in enumerate(ecEntropies):
            featEntropyD[(sName, ampN, "ecDNA_" + str(ind+1))] = etup

        if bfbStat:
            bfb_totalEnt, bfb_decompEnt, bfb_nEnt = decompositionComplexity(graphFile, cycleList, cycleCNs, segSeqD,
                                                                            bfb_cycle_inds, set())
            featEntropyD[(sName, ampN, "BFB_1")] = (bfb_totalEnt, bfb_decompEnt, bfb_nEnt)

        # write genes
        if args.report_genes:
            feat_genes = get_genes.extract_gene_list(sName, ampN, gene_lookup, args.report_genes, cycleList, segSeqD,
                                                     bfb_cycle_inds, ecIndexClusters, invalidInds, bfbStat, ecStat, ampClass)

            ftgd_list.append([sName, ampN, feat_genes])

        # store this additional information
        AMP_classifications.append((ampClass, ecStat, bfbStat, ecAmpliconCount))
        dvalues = [AMP_dvaluesDict[x] for x in categories]
        AMP_dvaluesList.append(dvalues)

        # edge classification
        edgeTypeCountD = defaultdict(float)
        if graphFile:
            posCycleLookup = buildPosCycleLookup(cycleList, segSeqD)
            bps = parseBPG(graphFile, args.add_chr_tag)
            for bp in bps:
                lCycles, rCycles = bpgEdgeToCycles(bp, posCycleLookup)
                # indices of left and right cycles on the discordant edges, and the index-ordered list of types
                resD = classifyConnections(lCycles, rCycles, cycleTypes)
                for k, v in resD.items():
                    edgeTypeCountD[mixLookups[k]] += v

            # norm the values
            eTCDSum = float(sum(edgeTypeCountD.values()))
            for k, v in edgeTypeCountD.items():
                edgeTypeCountD[k] = v / eTCDSum

        dvalues = [edgeTypeCountD[x] for x in mixing_cats]
        EDGE_dvaluesList.append(dvalues)

    # PLOTTING
    textCategories = ["No amp/Invalid", "Linear\namplification", "Trivial\ncycle", "Complex\nnon-cyclic",
                      "Complex\ncyclic", "BFB\nfoldback"]
    if args.plotStyle == "grouped":
        from radar_plotting import *

        print("plotting")
        make_classification_radar(textCategories, AMP_dvaluesList, args.o + "_amp_class", sampNames)
        make_classification_radar(mixing_cats, EDGE_dvaluesList, args.o + "_edge_class", sampNames)

    elif args.plotStyle == "individual":
        from radar_plotting import *

        print("plotting")
        for a, e, s in zip(AMP_dvaluesList, EDGE_dvaluesList, sampNames):
            print(textCategories, a)
            make_classification_radar(textCategories, [a[:len(textCategories)], ], args.o + "_" + s + "_amp_class",
                                      sampNames)
            make_classification_radar(mixing_cats, [e, ], args.o + "_" + s + "_edge_class", sampNames)

    #OUTPUT FILE WRITING
    print("writing output files")
    # Genes
    if args.report_genes:
        gene_extraction_outname = args.o + "_gene_list.tsv"
        get_genes.write_results(gene_extraction_outname, ftgd_list)

    # Feature entropy
    if args.report_complexity:
        with open(args.o + "_feature_entropy.tsv",'w') as outfile:
            outfile.write("sample\tamplicon\tfeature\ttotal_feature_entropy\tdecomp_entropy\tAmp_nseg_entropy\n")
            for k, vt in featEntropyD.items():
                ol = map(str, k + vt)
                outfile.write("\t".join(ol) + "\n")

    # Amplicon profiles
    with open(args.o + "_amplicon_classification_profiles.tsv", 'w') as outfile:
        oh = ["sample_name", "amplicon_number", "amplicon_decomposition_class", "ecDNA+", "BFB+", "ecDNA_amplicons"]
        if args.verbose_classification:
            oh += categories

        outfile.write("\t".join(oh) + "\n")
        for ind, sname in enumerate(sampNames):
            ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
            ampClass, ecStat, bfbStat, ecAmpliconCount = AMP_classifications[ind]
            ecOut = "Positive" if ecStat else "None detected"
            bfbOut = "Positive" if bfbStat else "None detected"
            ov = [sname.rsplit("_amplicon")[0], ampN, ampClass, ecOut, bfbOut, str(ecAmpliconCount)]
            if args.verbose_classification:
                ov+=[str(x) for x in AMP_dvaluesList[ind]]

            outfile.write("\t".join(ov) + "\n")


    # Edge profiles
    if args.verbose_classification:
        with open(args.o + "_edge_classification_profiles.tsv", 'w') as outfile:
            outfile.write("\t".join(["sample_name", "amplicon_number"] + mixing_cats) + "\n")
            for ind, sname in enumerate(sampNames):
                ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
                outfile.write(
                    "\t".join([sname.rsplit("_amplicon")[0], ampN] + [str(x) for x in EDGE_dvaluesList[ind]]) + "\n")

    print("done")
