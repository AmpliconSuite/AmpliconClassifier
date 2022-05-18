#!/usr/bin/env python3

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"
from amplicon_similarity import *
cn_cut = 4.5

# union of intervaltrees for each chromosome
def join_ivalds(ivald1, ivald2):
    joined_ivald = defaultdict(IntervalTree)
    ks1, ks2 = set(ivald1.keys()), set(ivald2.keys())
    for c in ks1.union(ks2):
        joined_ivald[c] = ivald1[c].union(ivald2[c])

    return joined_ivald


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute similarity between overlapping amplicon features, e.g. "
                                                 "ecDNAs, BFBs, etc.")
    parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, or GRCh38",
                        choices=["hg19", "GRCh37", "hg38", "GRCh38", "mm10", "GRCm38"], required=True)
    # parser.add_argument("-i", "--input", help="Path to list of AA amplicon files to use. Each line formatted as: \
    # samplename /path/to/sample_amplicon1_cycles.txt /path/to/sample_amplicon1_graph.txt", required=True)
    parser.add_argument("-f", "--feature_input", help="Path to list of AC feature bed files and corresponding graphs. "
                        "Input file for -f generated automatically by AC, 'feature_to_graph.txt'.", required=True)
    parser.add_argument("-o", help="Output filename prefix")
    parser.add_argument("--min_cn", type=float, help="Minimum CN to consider as amplification (default 4.5",
                        default=4.5)
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files",
                        action='store_true', default=False)
    parser.add_argument("--no_LC_filter", help="Do not filter low-complexity cycles. Not recommended to set this flag.",
                        action='store_true', default=False)
    parser.add_argument("--min_de", type=int, help="Set the minimum number of discordant edges in the amplicon "
                                                   "required to be considered for similarity (default 1).", default=1)
    parser.add_argument("--include_path_in_feature_name", help="Include path of file when reporting feature name "
                        "(useful for comparing against runs with similar names", action='store_true', default=False)
    parser.add_argument("--required_classifications", help="Which features to consider in the similarity calculation.",
                        choices=["ecDNA", "BFB", "Complex non-cyclic", "Linear amplification", "any"], nargs='+',
                        default="any")

    args = parser.parse_args()

    try:
        AC_SRC = os.environ["AC_SRC"]

    except KeyError:
        sys.stderr.write("$AC_SRC not found. Please first set AC_SRC bash variable using AC installation instructions "
                         "in Readme.\n")
        sys.exit(1)

    if args.ref == "hg38": args.ref = "GRCh38"
    add_chr_tag = args.add_chr_tag
    cn_cut = args.min_cn

    if not args.o:
        args.o = os.path.basename(args.input).rsplit(".")[0]

    outpre = args.o + "_feature_similarity/"
    if not os.path.exists(outpre):
        os.makedirs(outpre)

    try:
        AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + args.ref + "/"
    except KeyError:
        sys.stderr.write("$AA_DATA_REPO not set. Please see AA installation instructions.\n")
        sys.exit(1)

    lcD, cg5D = set_lcd(AA_DATA_REPO, args.no_LC_filter)

    # read bed files for each feature
    # amp2bed = defaultdict(lambda: defaultdict(str))
    feat_to_graph = {}
    feat_to_bed = {}
    total_feats = 0
    with open(args.feature_input) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            f, graphfile = fields[0], fields[1]
            featbasename = os.path.basename(f)
            featfields = featbasename.rsplit("_")
            feat = featfields[-3]
            if args.required_classifications == "any" or feat in args.required_classifications:
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

    with open(args.o + "_feature_similarity_scores.tsv", 'w') as outfile:
        outfile.write("Amp1\tAmp2\tSimilarityScore\tSimScorePercentile\tSimScorePvalue\tAsymmetricScore1\t"
                      "AsymmetricScore2\tGenomicSegmentScore1\tGenomicSegmentScore2\tBreakpointScore1\t"
                      "BreakpointScore2\tJaccardGenomicSegment\tJaccardBreakpoint\tNumSharedBPs\tAmp1NumBPs\t"
                      "Amp2NumBPs\tAmpOverlapLen\tAmp1AmpLen\tAmp2AmpLen\n")
        # for each pair, make the union bed
        outdata = []
        for x in pairs:
            s2a_graph = {}
            # rel_regions = join_ivalds(feat_to_bed[x[0]][1], feat_to_bed[x[1]][1])
            # read the two graph files and classify.
            graph0 = parseBPG(feat_to_graph[x[0]], feat_to_bed[x[0]][1], cn_cut, add_chr_tag, lcD, cg5D)
            graph1 = parseBPG(feat_to_graph[x[1]], feat_to_bed[x[1]][1], cn_cut, add_chr_tag, lcD, cg5D)
            if not graph0[0] or not graph0[1] or not graph1[0] or not graph1[1]:
                continue

            s2a_graph[x[0]] = graph0
            s2a_graph[x[1]] = graph1
            # for y,z in s2a_graph.items():
            #     print(y,z)

            compute_similarity(outfile, s2a_graph, [x], outdata)

        outdata.sort(key=lambda x: (x[2], x[1], x[0]), reverse=True)
        for l in outdata:
            outfile.write("\t".join([str(x) for x in l]) + "\n")