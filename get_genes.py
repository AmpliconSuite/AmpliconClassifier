#!/usr/bin/env python3

# This is imported by amplicon_classifier.py to get genes

from collections import defaultdict
import os

from intervaltree import IntervalTree


def parse_genes(gene_file):
    print("reading " + gene_file)
    t = defaultdict(IntervalTree)
    seenNames = set()
    with open(gene_file) as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split()
            if not fields:
                continue

            chrom, s, e, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            # parse the line and get the name
            propFields = {x.split("=")[0]: x.split("=")[1] for x in fields[-1].rstrip(";").split(";")}
            gname = propFields["Name"]
            is_other_feature = (gname.startswith("LOC") or gname.startswith("LINC") or gname.startswith("MIR"))
            if gname not in seenNames and not is_other_feature:
                seenNames.add(gname)
                t[chrom][s:e] = (gname, strand)

    print("read " + str(len(seenNames)) + " genes\n")
    return t


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


# take a list of 'feat_to_genes' dicts
def write_results(outname, ftg_list):
    with open(outname, 'w') as outfile:
        head = ["sample_name", "amplicon_number", "feature", "gene", "truncated"]
        outfile.write("\t".join(head) + "\n")
        for sname, anum, ftgd in ftg_list:
            for feat_name in sorted(ftgd.keys()):
                for gname in sorted(ftgd[feat_name].keys()):
                    truncs = [x for x in ["5p", "3p"] if x not in ftgd[feat_name][gname]]
                    ts = "_".join(truncs) if truncs else "None"
                    outfile.write("\t".join([sname, anum, feat_name, gname, ts]) + "\n")


#print all the intervals to bed files
def write_interval_beds(sname, ampN, feature_dict):
    outdir = "classification_bed_files/"
    os.makedirs(outdir,exist_ok=True)
    trim_sname = sname.rsplit("/")[-1]
    for feat_name, curr_fd in feature_dict.items():
        with open(outdir + trim_sname + "_" + ampN + "_" + feat_name + "_intervals.bed", 'w') as outfile:
            for chrom, ilist in curr_fd.items():
                if not chrom:
                    continue

                for i in ilist:
                    l = map(str, [chrom, i[0], i[1]])
                    outfile.write("\t".join(l) + "\n")


def extract_gene_list(sname, ampN, gene_lookup, classes_to_get, cycleList, segSeqD, bfb_cycle_inds, ecIndexClusters,
                      invalidInds, bfbStat, ecStat):

    feature_dict = {}
    if ("bfb" in classes_to_get or "both" in classes_to_get) and bfbStat:
        # collect unmerged genomic intervals comprising the feature
        bfb_interval_dict = defaultdict(list)
        for b_ind in bfb_cycle_inds:
            if b_ind not in invalidInds:
                for c_id in cycleList[b_ind]:
                    chrom, l, r = segSeqD[abs(c_id)]
                    bfb_interval_dict[chrom].append((l, r))

        feature_dict["BFB_1"] = bfb_interval_dict

    if ("ecdna" in classes_to_get or "both" in classes_to_get) and ecStat:
        # collect unmerged genomic intervals comprising the feature
        for amp_ind, ec_cycle_inds in enumerate(ecIndexClusters):
            ec_interval_dict = defaultdict(list)
            for e_ind in ec_cycle_inds:
                if e_ind not in invalidInds:
                    for c_id in cycleList[e_ind]:
                        chrom, l, r = segSeqD[abs(c_id)]
                        ec_interval_dict[chrom].append((l, r))

            feature_dict["ecDNA_" + str(amp_ind + 1)] = ec_interval_dict

    # merge all the intervals in each list of intervals
    tot_init_intervals = sum([len(ilist) for fd in feature_dict.values() for ilist in fd.values()])
    merge_intervals(feature_dict)
    tot_final_intervals = sum([len(ilist) for fd in feature_dict.values() for ilist in fd.values()])
    # print("Feature extraction: started with " + str(tot_init_intervals) + " unmerged intervals, finished with " + str(
    #     tot_final_intervals) + " intervals")

    write_interval_beds(sname, ampN, feature_dict)
    return get_genes_from_intervals(gene_lookup, feature_dict)
