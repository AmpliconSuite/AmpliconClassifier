from ac_util import *


# write a summary of the breakpoints
def summarize_breakpoints(graphf, add_chr_tag, feature_dict, lcD):
    linelist = [["chrom1", "pos1", "chrom2", "pos2", "sv_type", "read_support", "features", "orientation",
                 "pos1_flanking_coordinate", "pos2_flanking_coordinate"]]
    bps = bpg_edges(graphf, add_chr_tag, lcD)
    for bp in bps:
        c1, c2 = bp.lchrom, bp.rchrom
        p1, p2 = bp.lpos, bp.rpos
        ldir, rdir = bp.ldir, bp.rdir
        if c1 == c2:
            if p1 > p2:
                p1, p2 = p2, p1
                ldir, rdir = rdir, ldir

        if ldir == "+":
            p1_1before = p1 - 1
        else:
            p1_1before = p1 + 1

        if rdir == "+":
            p2_1before = p2 - 1
        else:
            p2_1before = p2 + 1

        if c1 != c2:
            etype = "interchromosomal"
        else:
            if ldir == "+" and rdir == "-":
                etype = "deletion-like"

            elif ldir == "-" and rdir == "+":
                etype = "duplication-like"

            else:
                if p2 - p1 < 25000:
                    etype = "foldback"
                else:
                    etype = "inversion"

        fl = []
        for f, intd in feature_dict.items():
            hitl, hitr = False, False
            for p in intd[c1]:
                if p[0] <= p1 <= p[1]:
                    hitl = True
                    break

            for p in intd[c2]:
                if p[0] <= p2 <= p[1]:
                    hitr = True
                    break

            if hitl and hitr:
                fl.append(f.rsplit("_")[0])

        if not fl:
            fl.append("None")

        fs = "|".join(fl)
        linelist.append([c1, p1, c2, p2, etype, bp.support, fs, ldir+rdir, p1_1before, p2_1before])

    return linelist


def get_gene_ends(gs, ge, strand, a, b):
    has5p, has3p = False, False
    if strand == "-":
        gs, ge = ge, gs

    if a <= gs <= b:
        has5p = True

    if a <= ge < b:
        has3p = True

    return has5p, has3p


def get_genes_from_intervals(gene_lookup, feature_dict, gseg_cn_d):
    ongene = os.path.dirname(os.path.realpath(__file__)) + "/resources/combined_oncogene_list.txt"
    ongene_set = set()
    with open(ongene) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            ongene_set.add(fields[0])

    feat_to_gene_trunc = defaultdict(lambda: defaultdict(set))
    feat_to_gene_cn = defaultdict(lambda: defaultdict(float))
    feat_to_ongene = defaultdict(lambda: defaultdict(bool))
    for feat_name, curr_fd in feature_dict.items():
        for chrom, intlist in curr_fd.items():
            for a, b in intlist:
                ogenes_ints = gene_lookup[chrom][a:b]
                for gp in ogenes_ints:
                    gname, strand = gp.data
                    has5p, has3p = get_gene_ends(gp.begin, gp.end, strand, a, b)
                    gsegs_hit = gseg_cn_d[chrom][gp.begin:gp.end]
                    gene_valid_cns = [x.data for x in gsegs_hit if x.end - x.begin > 1000]
                    if gene_valid_cns:
                        gene_cn = max(gene_valid_cns)
                    else:
                        gene_cn = "unknown"
                    feat_to_gene_cn[feat_name][gname] = gene_cn
                    if has5p:
                        feat_to_gene_trunc[feat_name][gname].add("5p")
                    if has3p:
                        feat_to_gene_trunc[feat_name][gname].add("3p")

                    if gname in ongene_set:
                        feat_to_ongene[feat_name][gname] = True

    return feat_to_gene_trunc, feat_to_gene_cn, feat_to_ongene


def get_gseg_cns(graphf, add_chr_tag):
    gseg_cn_d = defaultdict(IntervalTree)
    with open(graphf) as infile:
        for line in infile:
            fields = line.rsplit()
            if line.startswith("sequence"):
                c, s, e = fields[1].rsplit(":")[0], int(fields[1].rsplit(":")[1][:-1]), int(fields[2].rsplit(":")[1][:-1])+1
                cn = float(fields[3])
                if add_chr_tag and not c.startswith('chr'):
                    c = "chr" + c

                gseg_cn_d[c].addi(s, e, cn)

    return gseg_cn_d


def amplicon_len_and_cn(feature_dict, gseg_cn_d):
    prop_dict = {}
    for feat_name, curr_fd in feature_dict.items():
        if "unknown_" in feat_name:
            continue

        feat_size = 0
        cn_list = []
        for chrom, intlist in curr_fd.items():
            for a, b in intlist:
                feat_size+=(b-a)
                cns = gseg_cn_d[chrom][a:b]
                for d in cns:
                    if d.end - d.begin > 1000:
                        cn_list.append((d.data, min(d.end, b) - max(d.begin, a)))

        if not cn_list:
            continue

        cn_list.sort()
        tl = sum(x[1] for x in cn_list)
        mid = tl/2.0
        rs = 0
        median_cn = 0
        for cn, s in cn_list:
            rs+=s
            if rs >= mid:
                median_cn = cn
                break

        # max cn is the last element of the cn_list
        prop_dict[feat_name] = (feat_size, median_cn, cn_list[-1][0])

    return prop_dict


def amplicon_annotation(cycleList, segSeqD, bfb_cycle_inds, ecIndexClusters, invalidInds, bfbStat, ecStat, ampClass,
                        graphf, add_chr_tag, lcD):
    feature_dict = {}
    gseg_cn_d = get_gseg_cns(graphf, add_chr_tag)
    invalidSet = set(invalidInds)
    all_used = invalidSet.union(bfb_cycle_inds)
    used_segs = defaultdict(IntervalTree)
    graph_cns = get_graph_cns(graphf, add_chr_tag)
    other_class_c_inds = []
    if bfbStat:
        # collect unmerged genomic intervals comprising the feature
        bfb_interval_dict = defaultdict(list)
        for b_ind in bfb_cycle_inds:
            if b_ind not in invalidSet:
                for c_id in cycleList[b_ind]:
                    chrom, l, r = segSeqD[abs(c_id)]
                    if chrom:
                        bfb_interval_dict[chrom].append((l, r))
                        used_segs[chrom].addi(l, r+1)

        feature_dict["BFB_1"] = bfb_interval_dict

    if ecStat:
        # collect unmerged genomic intervals comprising the feature
        for amp_ind, ec_cycle_inds in enumerate(ecIndexClusters):
            ec_interval_dict = defaultdict(list)
            for e_ind in ec_cycle_inds:
                all_used.add(e_ind)
                if e_ind not in invalidSet:
                    for c_id in cycleList[e_ind]:
                        chrom, l, r = segSeqD[abs(c_id)]
                        if chrom:
                            used_segs[chrom].addi(l, r+1)
                            # chop out low cn regions
                            seg_t = IntervalTree([Interval(l, r+1)])
                            olapping_low_cns = [x for x in graph_cns[chrom][l:r+1] if x.data < 4]
                            for x in olapping_low_cns:
                                seg_t.chop(x.begin, x.end+1)

                            for x in seg_t:
                                ec_interval_dict[chrom].append((x.begin, x.end))

            feature_dict["ecDNA_" + str(amp_ind + 1)] = ec_interval_dict

    if ampClass != "No amp/Invalid":
        other_interval_dict = defaultdict(list)
        for o_ind in range(len(cycleList)):
            if o_ind not in all_used:
                for c_id in cycleList[o_ind]:
                    if abs(c_id) not in used_segs:
                        chrom, l, r = segSeqD[abs(c_id)]
                        if not used_segs[chrom][l:r] and not chrom is None:
                            other_interval_dict[chrom].append((l, r))
                            other_class_c_inds.append(o_ind)

        if not ecStat and not bfbStat:
            feature_dict[ampClass + "_1"] = other_interval_dict
        else:
            feature_dict["unknown_1"] = other_interval_dict

    merge_intervals(feature_dict)
    bpg_linelist = summarize_breakpoints(graphf, add_chr_tag, feature_dict, lcD)
    prop_dict = amplicon_len_and_cn(feature_dict, gseg_cn_d)
    return bpg_linelist, gseg_cn_d, other_class_c_inds, feature_dict, prop_dict
