#!/usr/bin/env python3

from ac_util import *


# This file is imported by amplicon_classifier.py to get genes

def merge_intervals(feature_dict):
    for item, usort_intd in feature_dict.items():
        for chrom, usort_ints in usort_intd.items():
            # sort ints
            sort_ints = sorted(usort_ints)
            # merge sorted ints
            mi = [sort_ints[0]]
            for ival in sort_ints[1:]:
                if ival[0] <= mi[-1][1]:
                    ui = (mi[-1][0], ival[1])
                    mi[-1] = ui

                else:
                    mi.append(ival)

            feature_dict[item][chrom] = mi


def get_gene_ends(gs, ge, strand, a, b):
    has5p, has3p = False, False
    if strand == "-":
        gs, ge = ge, gs

    if a <= gs <= b:
        has5p = True

    if a <= ge < b:
        has3p = True

    return has5p, has3p


def get_genes_from_intervals(gene_lookup, feature_dict):
    feat_to_genes = defaultdict(lambda: defaultdict(set))
    for feat_name, curr_fd in feature_dict.items():
        for chrom, intlist in curr_fd.items():
            for a, b in intlist:
                ogenes_ints = gene_lookup[chrom][a:b]
                for gp in ogenes_ints:
                    gname, strand = gp.data
                    has5p, has3p = get_gene_ends(gp.begin, gp.end, strand, a, b)
                    if has5p:
                        feat_to_genes[feat_name][gname].add("5p")
                    if has3p:
                        feat_to_genes[feat_name][gname].add("3p")

    return feat_to_genes


def extract_gene_list(sname, ampN, gene_lookup, classes_to_get, cycleList, segSeqD, bfb_cycle_inds, ecIndexClusters,
                      invalidInds, bfbStat, ecStat, ampClass):
    feature_dict = {}
    invalidSet = set(invalidInds)
    all_used = invalidSet.union(bfb_cycle_inds)
    if ("bfb" in classes_to_get or "all" in classes_to_get) and bfbStat:
        # collect unmerged genomic intervals comprising the feature
        bfb_interval_dict = defaultdict(list)
        for b_ind in bfb_cycle_inds:
            if b_ind not in invalidSet:
                for c_id in cycleList[b_ind]:
                    chrom, l, r = segSeqD[abs(c_id)]
                    bfb_interval_dict[chrom].append((l, r))

        feature_dict["BFB_1"] = bfb_interval_dict

    if ("ecdna" in classes_to_get or "all" in classes_to_get) and ecStat:
        # collect unmerged genomic intervals comprising the feature
        for amp_ind, ec_cycle_inds in enumerate(ecIndexClusters):
            ec_interval_dict = defaultdict(list)
            for e_ind in ec_cycle_inds:
                all_used.add(e_ind)
                if e_ind not in invalidSet:
                    for c_id in cycleList[e_ind]:
                        chrom, l, r = segSeqD[abs(c_id)]
                        ec_interval_dict[chrom].append((l, r))

            feature_dict["ecDNA_" + str(amp_ind + 1)] = ec_interval_dict

    if ("other" in classes_to_get or "all" in classes_to_get) and ampClass != "No amp/Invalid" and not ecStat and not bfbStat:
        other_interval_dict = defaultdict(list)
        for o_ind in range(len(cycleList)):
            if o_ind not in all_used:
                for c_id in cycleList[o_ind]:
                    chrom, l, r = segSeqD[abs(c_id)]
                    other_interval_dict[chrom].append((l, r))

        feature_dict["other_1"] = other_interval_dict

    # merge all the intervals in each list of intervals
    tot_init_intervals = sum([len(ilist) for fd in feature_dict.values() for ilist in fd.values()])
    merge_intervals(feature_dict)
    tot_final_intervals = sum([len(ilist) for fd in feature_dict.values() for ilist in fd.values()])
    # print("Feature extraction: started with " + str(tot_init_intervals) + " unmerged intervals, finished with " + str(
    #     tot_final_intervals) + " intervals")

    write_interval_beds(sname, ampN, feature_dict)
    return get_genes_from_intervals(gene_lookup, feature_dict)
