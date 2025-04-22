import logging

from intervaltree import IntervalTree

from ampclasslib.ac_util import *
from ampclasslib.ecDNA_context import *


def setup_logger(output_prefix):
    """Set up logging configuration"""
    log_file = "{}.log".format(output_prefix)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Setup file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)

    # Setup console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)

    # Setup logger
    logger = logging.getLogger('AmpliconClassifier')
    logger.setLevel(logging.INFO)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


# getting genes from the genes .gff file
def parse_genes(gene_file):
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

            try:
                ncbi_id = propFields["Accession"]
            except KeyError:
                ncbi_id = propFields["ID"]

            is_other_feature = (gname.startswith("LOC") or gname.startswith("LINC") or gname.startswith("MIR"))
            if gname not in seenNames and not is_other_feature:
                seenNames.add(gname)
                t[chrom][s:e] = (gname, strand, ncbi_id)

    return t


def bed_to_interval_dict(bedf, add_chr_tag):
    ivald = defaultdict(IntervalTree)
    with open(bedf) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if not fields:
                continue
            if add_chr_tag and not fields[0].startswith('chr'):
                fields[0] = 'chr' + fields[0]

            ivald[fields[0]].addi(int(fields[1]), int(fields[2]))

    return ivald


def parse_bpg(bpgf, add_chr_tag, lcD=None, subset_ivald=None, cn_cut=0, cg5D=None, min_de=0, min_de_size=0):
    """
    Parse breakpoint data from AA-formatted breakpoint graph files.

    Args:
        bpgf: Input breakpoint graph file path
        add_chr_tag: Whether to add 'chr' prefix to chromosome names
        lcD: Dictionary for chromosome positions of low-complexity elements
        subset_ivald: Interval dictionary for subsetting (default: None)
        cn_cut: Copy number cutoff (default: 0)
        cg5D: Dictionary for chromosome positioons of conserved gain regions (default: None)
        min_de: Minimum number of discordant edges in order to return parsed graph (default: 1)
        min_de_size: Minimum size for deletion events (default: 2500)

    Returns:
        Tuple of (list of breakpoints, segment tree)
    """

    if lcD is None:
        lcD = defaultdict(IntervalTree)

    if cg5D is None:
        cg5D = defaultdict(IntervalTree)

    bps = []
    segTree = defaultdict(IntervalTree)

    # Setup subsetting logic
    no_subset = subset_ivald is None or len(subset_ivald) == 0
    cn_cut_0 = (cn_cut <= 0 and no_subset)
    keepAll = (no_subset or cn_cut_0)

    # Edge parsing setup
    hom_len_index = None
    hom_seq_index = None

    with open(bpgf) as infile:
        for line in infile:
            if line.startswith("BreakpointEdge:"):
                fields = line.rstrip().rsplit()
                for ind, x in enumerate(fields):
                    if x.startswith("HomologySizeIfAvailable"):
                        hom_len_index = ind
                    elif x.startswith("Homology/InsertionSequence"):
                        hom_seq_index = ind

            elif line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                l, r = fields[1].rsplit("->")

                lchrom, lpos = l.rsplit(":")
                rchrom, rpos = r.rsplit(":")
                lpos, ldir = lpos[:-1], lpos[-1]
                rpos, rdir = rpos[:-1], rpos[-1]
                lpos, rpos = int(lpos), int(rpos)

                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                if lcD[lchrom][lpos] or lcD[rchrom][rpos]:
                    continue

                if not keepAll and not (subset_ivald[lchrom].overlaps(lpos) and subset_ivald[rchrom].overlaps(rpos)):
                    continue
                elif keepAll and not (segTree[lchrom].overlaps(lpos) or segTree[rchrom].overlaps(rpos)):
                    continue
                elif lchrom == rchrom and abs(lpos - rpos) < min_de_size and ldir != rdir:
                    continue

                cn = float(fields[2])
                support = int(fields[3]) if len(fields) > 3 else None

                # Get homology info if available
                homlen = homseq = "None"
                if not any([x is None for x in [hom_seq_index, hom_len_index]]):
                    homlen = fields[hom_len_index]
                    if homlen != "0":
                        homseq = fields[hom_seq_index]

                currBP = breakpoint(lchrom, lpos, rchrom, rpos, cn=cn, support=support,
                                    ldir=ldir, rdir=rdir, homlen=homlen, homseq=homseq)
                bps.append(currBP)

            elif line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                lchrom, lpos = fields[1].rsplit(":")
                lpos = int(lpos[:-1])
                rchrom, rpos = fields[2].rsplit(":")
                rpos = int(rpos[:-1]) + 1

                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = 'chr' + lchrom
                    rchrom = 'chr' + rchrom

                if lcD[lchrom][lpos:rpos] or cg5D[lchrom][lpos:rpos]:
                    continue

                if not keepAll and not subset_ivald[lchrom].overlaps(lpos, rpos):
                    continue

                cn = float(fields[3])
                if cn > cn_cut:
                    segTree[lchrom].addi(lpos, rpos, cn)

    if not bps and min_de > 0:
        return bps, defaultdict(IntervalTree)

    return bps, segTree


# write the number of ecDNA per sample
def write_ec_per_sample(outname, samp_to_ec_count, summary_map):
    # read the summary map and figure out what the sample names are
    all_samps = set()
    if isinstance(summary_map, str):
        with open(summary_map) as infile:
            for line in infile:
                all_samps.add(line.rsplit()[0])
    elif isinstance(summary_map, set):
        all_samps = summary_map  # assume it's a set

    seen_samps = set()
    with open(outname, 'w') as outfile:
        outfile.write("#sample\tecDNA_count\n")
        for s, c in samp_to_ec_count.items():
            outfile.write("\t".join([s, str(c)]) + "\n")
            seen_samps.add(s)

        unseen_samps = all_samps - seen_samps
        for s in sorted(unseen_samps):
            outfile.write("\t".join([s, "0"]) + "\n")


# take a list of 'feat_to_genes' dicts
def write_gene_results(outname, ftg_list):
    with open(outname, 'w') as outfile:
        head = ["sample_name", "amplicon_number", "feature", "gene", "gene_cn", "truncated", "is_canonical_oncogene", "ncbi_id"]
        outfile.write("\t".join(head) + "\n")
        for sname, anum, ampg_d in ftg_list:
            for feat_name in sorted(ampg_d.keys()):
                for curr_ag in sorted(ampg_d[feat_name].values(), key=lambda x: x.name):
                    trunc_string = curr_ag.get_truncation_status()
                    outfile.write("\t".join([sname, anum, feat_name, curr_ag.name, str(curr_ag.gene_cn), trunc_string,
                                             str(curr_ag.is_oncogene), curr_ag.ncbi_id]) + "\n")


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


def create_context_table(prefix, sname, ampN, feature_dict, samp_amp_to_graph, add_chr_tag=False):
    curr_contexts = []
    bedfile_dir = prefix + "_classification_bed_files/"
    cycles_dir = prefix + "_annotated_cycles_files/"
    trim_sname = sname.rsplit("/")[-1].rsplit("_amplicon")[0]
    for feat_name, curr_fd in feature_dict.items():
        intervals_fname = bedfile_dir + trim_sname + "_" + ampN + "_" + feat_name + "_intervals.bed"
        if feat_name.startswith("ecDNA"):
            if not curr_fd:
                context = "Unknown"

            else:
                gfile = samp_amp_to_graph[sname + "_" + ampN]
                base_graph_name = os.path.basename(gfile).rsplit("_graph.txt", 1)[0]
                cfile = cycles_dir + base_graph_name + "_annotated_cycles.txt"
                context = fetch_context(gfile, cfile, bed_file=intervals_fname, add_chr_tag=add_chr_tag)

            curr_contexts.append((trim_sname + "_" + ampN + "_" + feat_name, context))

    return curr_contexts


# write a cycles file with the cycles -> some corrected
def write_annotated_corrected_cycles_file(prefix, outname, cycleList, cycleCNs, segSeqD, bfb_cycle_inds, ecIndexClusters,
                                          invalidInds, rearrCycleInds):
    offset_d = defaultdict(lambda: 1)
    offset_d[0] = 0
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
            clen = sum([segSeqD[abs(x)][2] - segSeqD[abs(x)][1] + offset_d[x] for x in cyc])
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
                  samp_amp_to_graph, prop_list, summary_map, logname="AmpliconClassifier"):

    logger = logging.getLogger(logname)
    # Genes
    gene_extraction_outname = args.o + "_gene_list.tsv"
    write_gene_results(gene_extraction_outname, ftgd_list)
    ecDNA_count_outname = args.o + "_ecDNA_counts.tsv"
    write_ec_per_sample(ecDNA_count_outname, samp_to_ec_count, summary_map)

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

    feat_type_counts = defaultdict(int)

    for ind, (sname, bpg_linelist, feature_dict, prop_dict) in enumerate(zip(sampNames, bpgi_list, fd_list, prop_list)):
        ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
        write_bpg_summary(args.o, sname, ampN, bpg_linelist)
        write_interval_beds(args.o, sname, ampN, feature_dict, samp_amp_to_graph, f2gf)
        write_basic_properties(feat_basic_propf, sname, ampN, prop_dict)
        for feat_type_and_ind in feature_dict.keys():
            feat_type = feat_type_and_ind.rsplit("_")[0]
            feat_type_counts[feat_type] += 1

    f2gf.close()

    # Check if there are any feature types first
    logger.info("\nFeature Type Counts:")
    if feat_type_counts:
        # Get max length of feature type names for nice alignment
        max_name_length = max(len(feat_type) for feat_type in feat_type_counts.keys())

        # Log header
        logger.info("-" * (max_name_length + 10))  # Line separator

        # Log each count with aligned formatting
        for feat_type, count in sorted(feat_type_counts.items()):
            if feat_type == "unknown":
                continue
            logger.info("{:<{width}} : {:>5}".format(feat_type, count, width=max_name_length))

        logger.info("-" * (max_name_length + 10))  # Line separator
    else:
        # Handle empty feature case
        logger.info("-" * 20)  # Default line separator
        logger.info("No focal amp features identified")
        logger.info("-" * 20)  # Default line separator

    # report ecDNA context
    context_filename = args.o + "_ecDNA_context_calls.tsv"
    contexts = []
    for ind, (sname, feature_dict) in enumerate(zip(sampNames, fd_list)):
        ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
        curr_contexts = create_context_table(args.o, sname, ampN, feature_dict, samp_amp_to_graph, add_chr_tag=args.add_chr_tag)
        contexts.extend(curr_contexts)

    with open(context_filename, 'w') as context_outfile:
        for x in contexts:
            context_outfile.write("\t".join(x) + "\n")

    # Edge profiles
    if args.verbose_classification:
        with open(args.o + "_edge_classification_profiles.tsv", 'w') as outfile:
            outfile.write("\t".join(["sample_name", "amplicon_number"] + mixing_cats) + "\n")
            for ind, sname in enumerate(sampNames):
                ampN = cyclesFiles[ind].rstrip("_cycles.txt").rsplit("_")[-1]
                outfile.write(
                    "\t".join([sname.rsplit("_amplicon")[0], ampN] + [str(x) for x in EDGE_dvaluesList[ind]]) + "\n")
