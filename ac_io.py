from ac_annotation import *


# getting genes
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
            try:
                gname = propFields["Name"]
            except KeyError:
                gname = propFields["gene_name"]
            is_other_feature = (gname.startswith("LOC") or gname.startswith("LINC") or gname.startswith("MIR"))
            if gname not in seenNames and not is_other_feature:
                seenNames.add(gname)
                t[chrom][s:e] = (gname, strand)

    print("read " + str(len(seenNames)) + " genes\n")
    return t


# write the number of ecDNA per sample
def write_ec_per_sample(outname, samp_to_ec_count):
    with open(outname, 'w') as outfile:
        outfile.write("#sample\tecDNA_count\n")
        for s, c in samp_to_ec_count.items():
            outfile.write("\t".join([s, str(c)]) + "\n")


# take a list of 'feat_to_genes' dicts
def write_gene_results(outname, ftg_list):
    with open(outname, 'w') as outfile:
        head = ["sample_name", "amplicon_number", "feature", "gene", "gene_cn", "truncated", "is_canonical_oncogene"]
        outfile.write("\t".join(head) + "\n")
        for sname, anum, truncd, cnd, ogd in ftg_list:
            for feat_name in sorted(truncd.keys()):
                for gname in sorted(truncd[feat_name].keys()):
                    truncs = [x for x in ["5p", "3p"] if x not in truncd[feat_name][gname]]
                    gene_cn = str(cnd[feat_name][gname])
                    inongene = str(ogd[feat_name][gname])
                    ts = "_".join(truncs) if truncs else "None"
                    outfile.write("\t".join([sname, anum, feat_name, gname, gene_cn, ts, inongene]) + "\n")


def write_basic_properties(feat_basic_propf, sname, ampN, prop_dict):
    for feat_name, props in prop_dict.items():
        borderline_flag = ""
        if (feat_name.startswith("ecDNA") or feat_name.startswith("BFB")) and props[-1] < 8:
            borderline_flag+="LowCN"
        elif props[-1] < 5:
            borderline_flag+="LowCN"

        if not borderline_flag:
            borderline_flag = "None"

        feat_basic_propf.write("\t".join([sname + "_" + ampN + "_" + feat_name, ] + [str(x) for x in props] + [borderline_flag,]) + "\n")


# write a summary of the breakpoints in the graph
def write_bpg_summary(prefix, sname, ampN, bpg_linelist):
    outdir = prefix + "_SV_summaries/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    trim_sname = sname.rsplit("/")[-1].rsplit("_amplicon")[0]
    full_fname = outdir + trim_sname + "_" + ampN + "_SV_summary.tsv"
    with open(full_fname, 'w') as outfile:
        for line in bpg_linelist:
            outfile.write("\t".join([str(x) for x in line]) + "\n")


# print all the intervals to bed files
def write_interval_beds(prefix, sname, ampN, feature_dict, samp_amp_to_graph, f2gf):
    # the matching and filename stuff is already done in amplicon_classifier. consider refactoring to make use of that
    outdir = prefix + "_classification_bed_files/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    trim_sname = sname.rsplit("/")[-1].rsplit("_amplicon")[0]
    for feat_name, curr_fd in feature_dict.items():
        if not curr_fd:
            continue

        full_fname = outdir + trim_sname + "_" + ampN + "_" + feat_name + "_intervals.bed"
        if not feat_name.startswith("unknown"):
            f2gf.write(full_fname + "\t" + samp_amp_to_graph[sname + "_" + ampN] + "\n")
        with open(full_fname, 'w') as outfile:
            for chrom, ilist in curr_fd.items():
                if not chrom:
                    continue

                for i in ilist:
                    l = map(str, [chrom, i[0], i[1]])
                    outfile.write("\t".join(l) + "\n")


# write a cycles file with the cycles -> some corrected
def write_annotated_corrected_cycles_file(prefix, outname, cycleList, cycleCNs, segSeqD, bfb_cycle_inds, ecIndexClusters,
                                          invalidInds, rearrCycleInds):
    outdir = prefix + "_annotated_cycles_files/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    with open(outdir + outname, 'w') as outfile:
        outfile.write("List of cycle segments\n")
        seg_ind = 0
        for k in sorted(segSeqD.keys()):
            if k != 0:
                seg_ind += 1
                ll = "\t".join(["Segment", str(seg_ind), segSeqD[k][0], str(segSeqD[k][1]), str(segSeqD[k][2])])
                outfile.write(ll + "\n")

        for ind, cyc in enumerate(cycleList):
            ccn = cycleCNs[ind]
            clen = sum([segSeqD[abs(x)][2] - segSeqD[abs(x)][1] for x in cyc])
            cl = ",".join([str(abs(x)) + "-" if x < 0 else str(abs(x)) + "+" for x in cyc])
            acclass = ""
            isCyclic = cyc[0] != 0
            for ec_i, ec_cycset in enumerate(ecIndexClusters):
                if ind in ec_cycset:
                    acclass += "ecDNA-like_"

            if ind in bfb_cycle_inds:
                acclass += "BFB-like_"

            if ind in invalidInds:
                acclass += "Invalid"

            acclass = acclass.rstrip("_")
            if not acclass:
                if ind in rearrCycleInds or isCyclic:
                    acclass = "Rearranged"
                else:
                    acclass = "Linear"

            l1 = ";".join(["Cycle=" + str(ind+1), "Copy_count=" + str(ccn), "Length=" + str(clen),
                           "IsCyclicPath=" + str(isCyclic), "CycleClass=" + acclass, "Segments=" + cl])
            outfile.write(l1 + "\n")


def write_outputs(args, ftgd_list, ftci_list, bpgi_list, featEntropyD, categories, sampNames, cyclesFiles,
                  AMP_classifications, AMP_dvaluesList, mixing_cats, EDGE_dvaluesList, samp_to_ec_count, fd_list,
                  samp_amp_to_graph, prop_list):
    # Genes
    gene_extraction_outname = args.o + "_gene_list.tsv"
    write_gene_results(gene_extraction_outname, ftgd_list)
    ecDNA_count_outname = args.o + "_ecDNA_counts.tsv"
    write_ec_per_sample(ecDNA_count_outname, samp_to_ec_count)

    # Feature entropy
    if args.report_complexity:
        with open(args.o + "_feature_entropy.tsv", 'w') as outfile:
            outfile.write("sample\tamplicon\tfeature\ttotal_feature_entropy\tdecomp_entropy\tAmp_nseg_entropy\n")
            sorted_keys = sorted(featEntropyD.keys())
            for k in sorted_keys:
                vt = featEntropyD[k]
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
                ov += [str(x) for x in AMP_dvaluesList[ind]]

            outfile.write("\t".join(ov) + "\n")

    for x in ftci_list:
        outname, cycleList, cycleCNs, segSeqD, bfb_cycle_inds, ecIndexClusters, invalidInds, rearrCycleInds = x
        write_annotated_corrected_cycles_file(args.o, outname, cycleList, cycleCNs, segSeqD, bfb_cycle_inds,
                                              ecIndexClusters, invalidInds, rearrCycleInds)

    feat_basic_propf = open(args.o + "_feature_basic_properties.tsv", 'w')
    prop_head = ["feature_ID", "captured_region_size_bp", "median_feature_CN", "max_feature_CN", "borderline_flag"]
    feat_basic_propf.write("\t".join(prop_head) + "\n")
    f2gf = open(args.o + "_features_to_graph.txt", 'w')

    for ind, (sname, bpg_linelist, feature_dict, prop_dict) in enumerate(zip(sampNames, bpgi_list, fd_list, prop_list)):
        ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
        write_bpg_summary(args.o, sname, ampN, bpg_linelist)
        write_interval_beds(args.o, sname, ampN, feature_dict, samp_amp_to_graph, f2gf)
        write_basic_properties(feat_basic_propf, sname, ampN, prop_dict)

    f2gf.close()

    # Edge profiles
    if args.verbose_classification:
        with open(args.o + "_edge_classification_profiles.tsv", 'w') as outfile:
            outfile.write("\t".join(["sample_name", "amplicon_number"] + mixing_cats) + "\n")
            for ind, sname in enumerate(sampNames):
                ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
                outfile.write(
                    "\t".join([sname.rsplit("_amplicon")[0], ampN] + [str(x) for x in EDGE_dvaluesList[ind]]) + "\n")
