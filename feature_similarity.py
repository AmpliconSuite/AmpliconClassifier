#!/usr/bin/env python3

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"

from ampclasslib.amplicon_similarity import *
from concurrent.futures import ProcessPoolExecutor, as_completed

cn_cut = 4.5
FEATURE_SIMILARITY_WORKER_CONTEXT = None


# union of intervaltrees for each chromosome
def join_ivalds(ivald1, ivald2):
    joined_ivald = defaultdict(IntervalTree)
    ks1, ks2 = set(ivald1.keys()), set(ivald2.keys())
    for c in ks1.union(ks2):
        joined_ivald[c] = ivald1[c].union(ivald2[c])

    return joined_ivald


def _init_feature_similarity_worker(feat_to_graph, feat_to_bed, add_chr_tag, lcD, cg5D, cn_cut_value, min_de_value):
    global FEATURE_SIMILARITY_WORKER_CONTEXT
    FEATURE_SIMILARITY_WORKER_CONTEXT = {
        "feat_to_graph": feat_to_graph,
        "feat_to_bed": feat_to_bed,
        "add_chr_tag": add_chr_tag,
        "lcD": lcD,
        "cg5D": cg5D,
        "cn_cut": cn_cut_value,
        "min_de": min_de_value,
    }


def _compute_feature_pair_similarity(pair):
    context = FEATURE_SIMILARITY_WORKER_CONTEXT
    feat_to_graph = context["feat_to_graph"]
    feat_to_bed = context["feat_to_bed"]
    feature_id_1, feature_id_2 = pair
    graph_path_1 = feat_to_graph[feature_id_1]
    graph_path_2 = feat_to_graph[feature_id_2]
    if graph_path_1 == graph_path_2:
        return []

    graph1 = parse_bpg(
        graph_path_1, context["add_chr_tag"], context["lcD"], subset_ivald=feat_to_bed[feature_id_1][1],
        cn_cut=context["cn_cut"], cg5D=context["cg5D"], min_de=context["min_de"], min_de_size=min_de_size
    )
    graph2 = parse_bpg(
        graph_path_2, context["add_chr_tag"], context["lcD"], subset_ivald=feat_to_bed[feature_id_2][1],
        cn_cut=context["cn_cut"], cg5D=context["cg5D"], min_de=context["min_de"], min_de_size=min_de_size
    )
    if (not graph1[0] and context["min_de"] > 0) or not graph1[1]:
        return []
    if (not graph2[0] and context["min_de"] > 0) or not graph2[1]:
        return []

    outdata = []
    compute_similarity({feature_id_1: graph1, feature_id_2: graph2}, [(feature_id_1, feature_id_2)], outdata)
    return outdata


def compute_pairwise_feature_similarity(feat_to_graph, feat_to_bed, pairs, add_chr_tag=False, lcD=None, cg5D=None,
                                        cn_cut_value=0, min_de_value=0, threads=1):
    if lcD is None:
        lcD = defaultdict(IntervalTree)

    worker_count = max(1, int(threads or 1))
    outdata = []
    if worker_count == 1 or len(pairs) <= 1:
        _init_feature_similarity_worker(feat_to_graph, feat_to_bed, add_chr_tag, lcD, cg5D, cn_cut_value,
                                        min_de_value)
        for pair in pairs:
            outdata.extend(_compute_feature_pair_similarity(pair))
    else:
        worker_count = min(worker_count, len(pairs))
        with ProcessPoolExecutor(
            max_workers=worker_count,
            initializer=_init_feature_similarity_worker,
            initargs=(feat_to_graph, feat_to_bed, add_chr_tag, lcD, cg5D, cn_cut_value, min_de_value)
        ) as executor:
            futures = [executor.submit(_compute_feature_pair_similarity, pair) for pair in pairs]
            for future in as_completed(futures):
                outdata.extend(future.result())

    outdata.sort(key=lambda x: (x[2], x[1], x[0]), reverse=True)
    return outdata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute similarity between overlapping amplicon features, e.g. "
                                                 "ecDNAs, BFBs, etc.")
    parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, or GRCh38",
                        choices=["hg19", "GRCh37", "hg38", "GRCh38", "mm10", "GRCm38", "GRCh38_viral"], required=True)
    parser.add_argument("-f", "--feature_input", help="Path to list of AC feature bed files and corresponding graphs. "
                        "Input file for -f generated automatically by AC, '[sample]_features_to_graph.txt'.",
                        required=True)
    parser.add_argument("-o", help="Output filename prefix")
    parser.add_argument("--min_cn", type=float, help="Minimum CN to consider as amplification from regions previously "
                        "assigned in feature bed file (default 0).", default=0)
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files",
                        action='store_true', default=False)
    parser.add_argument("--no_LC_filter", help="Do not filter low-complexity cycles. Not recommended to set this flag.",
                        action='store_true', default=False)
    parser.add_argument("--min_de", type=int, help="Set the minimum number of discordant edges in the amplicon "
                                                   "required to be considered for similarity (default 0).", default=0)
    parser.add_argument("--include_path_in_feature_name", help="Include path of file when reporting feature name "
                        "(useful for comparing against runs with similar names", action='store_true', default=False)
    parser.add_argument("--required_classifications", help="Which features to consider in the similarity calculation. "
                        "(default 'any')", choices=["ecDNA", "BFB", "CNC", "Linear", "any"], nargs='+', default="any")
    parser.add_argument("--threads", help="Number of worker processes for pairwise similarity scoring (default 1).",
                        type=int, default=1)

    args = parser.parse_args()

    if args.threads < 1:
        sys.stderr.write("--threads must be >= 1\n")
        sys.exit(1)

    if args.ref == "hg38": args.ref = "GRCh38"
    elif args.ref == "GRCm38": args.ref = "mm10"

    try:
        AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + args.ref + "/"
    except KeyError:
        sys.stderr.write("$AA_DATA_REPO not set. Please see AA installation instructions.\n")
        sys.exit(1)

    lcD, cg5D = set_lcd(AA_DATA_REPO, args.no_LC_filter)

    if "any" in args.required_classifications:
        required_classes = {"ecDNA", "BFB", "Complex-non-cyclic", "Linear"}
    else:
        required_classes = set(args.required_classifications)

    if "CNC" in required_classes:
        required_classes.remove("CNC")
        required_classes.add("Complex-non-cyclic")

    if "Linear" in required_classes:
        required_classes.remove("Linear")
        required_classes.add("Linear")

    print("Required classifications set to: " + str(required_classes))

    add_chr_tag = args.add_chr_tag
    cn_cut = args.min_cn

    if not args.o:
        args.o = os.path.basename(args.feature_input).rsplit(".")[0]
        if args.o.endswith("features_to_graph"):
            args.o = args.o.rsplit("_features_to_graph", 1)[0]

    odir = "/".join(args.o.rsplit("/")[:-1])
    if odir != "" and not os.path.exists(odir):
        os.makedirs(odir)

    # read bed files for each feature
    # amp2bed = defaultdict(lambda: defaultdict(str))
    feat_to_graph = {}
    feat_to_bed = {}
    total_feats = 0
    with open(args.feature_input) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            if not fields:
                continue
            f, graphfile = fields[0], fields[1]
            featbasename = os.path.basename(f)
            featfields = featbasename.rsplit("_")
            feat = featfields[-3]
            if feat in required_classes:
                # fullpath = classBedDir + "/" + f
                # amp2bed[amp][feat] = bed_to_interval_dict(fullpath)
                featname = f.rsplit("_intervals.bed")[0]
                if not args.include_path_in_feature_name:
                    featname = featname.rsplit("/")[-1]

                feat_to_graph[featname] = graphfile
                feat_to_bed[featname] = (None, bed_to_interval_dict(f, add_chr_tag))
                total_feats += 1

    print("Total of " + str(total_feats) + " features.")

    # get pairs
    pairs = get_pairs(feat_to_bed)
    print("Total of " + str(len(pairs)) + " pairs of features to compare.")
    if args.threads > 1 and len(pairs) > 1:
        print("Computing feature similarity with " + str(min(args.threads, len(pairs))) + " worker processes.")

    with open(args.o + "_feature_similarity_scores.tsv", 'w') as outfile:
        outfile.write("Amp1\tAmp2\tSimilarityScore\tSimScorePercentile\tSimScorePvalue\tAsymmetricScore1\t"
                      "AsymmetricScore2\tGenomicSegmentScore1\tGenomicSegmentScore2\tBreakpointScore1\t"
                      "BreakpointScore2\tJaccardGenomicSegment\tJaccardBreakpoint\tNumSharedBPs\tAmp1NumBPs\t"
                      "Amp2NumBPs\tAmpOverlapLen\tAmp1AmpLen\tAmp2AmpLen\n")
        outdata = compute_pairwise_feature_similarity(
            feat_to_graph, feat_to_bed, pairs, add_chr_tag=add_chr_tag, lcD=lcD, cg5D=cg5D,
            cn_cut_value=cn_cut, min_de_value=args.min_de, threads=args.threads
        )
        for l in outdata:
            outfile.write("\t".join([str(x) for x in l]) + "\n")
