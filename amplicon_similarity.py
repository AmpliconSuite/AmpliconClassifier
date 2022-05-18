#!/usr/bin/env python3

__author__ = "Jens Luebeck (jluebeck [at] ucsd.edu)"

import argparse
import bisect
from collections import defaultdict
import os
import sys

from intervaltree import IntervalTree
from scipy.stats import beta

add_chr_tag = False
bpweight = 0.75
d = 250
cn_cut = 4.5
min_de = 1
min_de_size = 2500
# negative log likelihood model fit similarity score model parameters from null distribution
alpha0, beta0 = 1.0407843068316383, 7.0222490956961359
'''
input:
tab-separated file listing for each amplicon/sample 
amplicon_name    amplicon_graph.txt    amplicon_cycles.txt


output:
pair wise similarity between all amplicons on a single line, for all lines in the input file

'''

class breakpoint(object):
    def __init__(self, lchrom, lpos, rchrom, rpos, cn):
        self.lchrom = lchrom
        self.lpos = lpos
        self.rchrom = rchrom
        self.rpos = rpos
        self.cn = cn

    def to_string(self):
        return self.lchrom + ":" + str(self.lpos) + " | " + self.rchrom + ":" + str(self.rpos) + "\t" + str(self.cn)

    def d_similar(self, bp2, d):
        bp2_chrom_set = {bp2.lchrom, bp2.rchrom}
        if self.lchrom not in bp2_chrom_set or self.rchrom not in bp2_chrom_set:
            return False

        sbp1 = sorted([(self.lchrom, self.lpos), (self.rchrom, self.rpos)])
        sbp2 = sorted([(bp2.lchrom, bp2.lpos), (bp2.rchrom, bp2.rpos)])

        if sbp1[0][0] == sbp2[0][0] and sbp1[1][0] == sbp2[1][0]:
            if abs(sbp1[0][1] - sbp2[0][1]) + abs(sbp1[1][1] - sbp2[1][1]) < d:
                return True

        return False


def amps_overlap(gdoi1, gdoi2):
    for c1, it1 in gdoi1.items():
        it2 = gdoi2[c1]
        if not it2:
            continue

        for int1 in it1:
            if int1.data is None or int1.data > cn_cut:
                for int2 in it2[int1.begin:int1.end]:
                    if int2.data is None or int2.data > cn_cut:
                        return True

    return False


def nucSimilarity(g1, g2):
    g1_bp = 0
    obp = 0
    for c1, it1 in g1.items():
        for t1 in it1:
            g1_bp += (t1.end - t1.begin)
            for t2 in g2[c1][t1.begin:t1.end]:
                obp += (min(t1.end, t2.end) - max(t1.begin, t2.begin))

    return obp/g1_bp


def bp_dist(bplist1, bplist2, d):
    intersectCount = 0.0
    for bp1 in bplist1:
        for bp2 in bplist2:
            if bp1.d_similar(bp2,d):
                intersectCount += 1
                break

    tE = max(len(bplist1), 1)
    return intersectCount/tE


def jaccard_sim_seq(g1, g2):
    g1_bp = 0
    g2_bp = 0
    obp = 0
    for c1, it1 in g1.items():
        for t1 in it1:
            g1_bp += (t1.end - t1.begin)
            for t2 in g2[c1][t1.begin:t1.end]:
                obp += (min(t1.end, t2.end) - max(t1.begin, t2.begin))

    for c2, it2 in g2.items():
        for t2 in it2:
             g2_bp += (t2.end - t2.begin)

    jsim = obp/(g1_bp + g2_bp - obp)
    return jsim, obp, g1_bp, g2_bp


def compute_num_shared_bps(bplist1, bplist2, d):
    if not bplist1 and not bplist1:
        return 0

    intersectCount = 0.0
    for bp1 in bplist1:
        for bp2 in bplist2:
            if bp1.d_similar(bp2, d):
                intersectCount += 1
                break

    return intersectCount


def jaccard_sim_bp(bplist1, bplist2, d):
    if not bplist1 and not bplist1:
        return 0, 0

    intersectCount = compute_num_shared_bps(bplist1, bplist2, d)
    jsim = intersectCount/(len(bplist1) + len(bplist2) - intersectCount)
    return jsim, intersectCount


def asymmetric_score(bplist1, bplist2, st1, st2, d):
    ns = nucSimilarity(st1, st2)
    bpd = bp_dist(bplist1, bplist2, d)
    S = (1 - bpweight) * ns + (bpweight) * bpd
    return S, ns, bpd


def symScore(bplist_a, bplist_b, st_a, st_b, d):
    as1, nS1, bpd1 = asymmetric_score(bplist_a, bplist_b, st_a, st_b, d)
    as2, nS2, bpd2 = asymmetric_score(bplist_b, bplist_a, st_b, st_a, d)
    return (as1 + as2)/2, [as1, as2, nS1, nS2, bpd1, bpd2]


def score_to_pval(s):
    return 1-beta.cdf(s, alpha0, beta0)


def score_to_percentile(s, background_scores):
    ind = bisect.bisect(background_scores, s)
    pval = ind*100.0/len(background_scores)
    return pval


def get_pairs(s2a_graph, required_amplicon_ID=None):
    pairs = []
    graph_data = []

    dat = list(s2a_graph.items())
    for i1 in range(len(dat)):
        s1, gpair1 = dat[i1]
        _, segTree = gpair1
        graph_data.append(segTree)
        for i2 in range(i1):
            s2, _ = dat[i2]
            if required_amplicon_ID and (not s1.endswith(required_amplicon_ID) and not s2.endswith(required_amplicon_ID)):
                continue

            if amps_overlap(graph_data[i1], graph_data[i2]):
                pairs.append((s1, s2))

    print(str(len(pairs)) + " overlapping pairs found.")
    return pairs


def buildLCDatabase(mappabilityFile):
    lcD = defaultdict(IntervalTree)
    with open(mappabilityFile) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            chrom, s, e = fields[0], int(fields[1]), int(fields[2])
            if e - s > 7500:
                lcD[chrom].addi(s, e)

    return lcD


def build_CG5_database(cg5_file):
    cg5D = defaultdict(IntervalTree)
    with open(cg5_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            chrom, s, e = fields[0], int(fields[1]), int(fields[2])
            cg5D[chrom].addi(s, e)

    return cg5D


def parseBPG(bpgf, subset_ivald, cn_cut, add_chr_tag, lcD, cg5D):
    bps = []
    keepAll = (len(subset_ivald) == 0 or cn_cut == 0)
    segTree = defaultdict(IntervalTree)
    with open(bpgf) as infile:
        for line in infile:
            # if line.startswith("discordant") or line.startswith("concordant"):
            if line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                l, r = fields[1].rsplit("->")
                lchrom, lpos = l[:-1].rsplit(":")
                rchrom, rpos = r[:-1].rsplit(":")
                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = "chr" + lchrom
                    rchrom = "chr" + rchrom

                lpos, rpos = int(lpos), int(rpos)
                if lcD[lchrom][lpos] or lcD[rchrom][rpos]:
                    continue

                elif not keepAll and not (subset_ivald[lchrom].overlaps(lpos) and subset_ivald[rchrom].overlaps(rpos)):
                    continue

                elif keepAll and not (segTree[lchrom].overlaps(lpos) or segTree[rchrom].overlaps(rpos)):
                    continue

                elif lchrom == rchrom and abs(lpos - rpos) < min_de_size:
                    continue

                cn = float(fields[2])
                currBP = breakpoint(lchrom, lpos, rchrom, rpos, cn)
                bps.append(currBP)

            elif line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                lchrom, lpos = fields[1].rsplit(":")
                lpos = int(lpos[:-1])
                rchrom, rpos = fields[2].rsplit(":")
                rpos = int(rpos[:-1])+1

                if add_chr_tag and not lchrom.startswith('chr'):
                    lchrom = 'chr' + lchrom
                    rchrom = 'chr' + rchrom

                if lcD[lchrom][lpos:rpos] or cg5D[lchrom][lpos:rpos]:
                    continue

                elif not keepAll and not subset_ivald[lchrom].overlaps(lpos, rpos):
                    continue

                cn = float(fields[3])
                if cn > cn_cut:
                    segTree[lchrom].addi(lpos, rpos, cn)

    if not bps:
        return bps, defaultdict(IntervalTree)

    return bps, segTree


def readFlist(filelist, region_subset_ivald, lcD, cg5D, cn_cut=0, fullname=False, required_classes=None, a2class=None):
    if required_classes is None:
        required_classes = set()
    if a2class is None:
        a2class = defaultdict(set)

    s_to_amp_graph_lookup = {}
    with open(filelist) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                if len(fields) > 1:
                    if fullname:
                        sname = fields[2].rsplit("_graph.txt")[0]
                    else:
                        sname = fields[2].rsplit("/")[-1].rsplit("_graph.txt")[0]
                        if sname.startswith("AA"):
                            sname = fields[0]

                    if not required_classes or (required_classes and a2class[sname].intersection(required_classes)):
                        s_to_amp_graph_lookup[sname] = parseBPG(fields[2], region_subset_ivald, cn_cut, add_chr_tag,
                                                                lcD, cg5D)

    return s_to_amp_graph_lookup


def read_classifications(classification_file):
    a2class = defaultdict(set)
    with open(classification_file) as infile:
        headfields = next(infile).rstrip().rsplit("\t")
        for line in infile:
            ld = dict(zip(headfields, line.rstrip().rsplit("\t")))
            ampid = ld["sample_name"] + "_" + ld["amplicon_number"]
            if ld["ecDNA+"] == "Positive":
                a2class[ampid].add("ecDNA")
            if ld["BFB+"] == "Positive":
                a2class[ampid].add("BFB")

            a2class[ampid].add(ld["amplicon_decomposition_class"])

    return a2class


def readBackgroundScores(bgsfile):
    with open(bgsfile) as infile:
        scores = [float(x) for x in next(infile).rstrip().rsplit("\t")]
        return scores


def bed_to_interval_dict(bedf, add_chr_tag):
    ivald = defaultdict(IntervalTree)
    with open(bedf) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if add_chr_tag and not fields[0].startswith('chr'):
                fields[0] = 'chr' + fields[0]

            ivald[fields[0]].addi(int(fields[1]), int(fields[2]))

    return ivald


def compute_similarity(outfile, s2a_graph, pairs, outdata):
    bgsfile = os.path.dirname(os.path.realpath(__file__)) + "/resources/sorted_background_scores.txt"
    background_scores = readBackgroundScores(bgsfile)
    for a, b in pairs:
        bplist_a, st_a = s2a_graph[a]
        bplist_b, st_b = s2a_graph[b]
        s, featList = symScore(bplist_a, bplist_b, st_a, st_b, d)
        jaccard_seq, amp_olap_len, amp_a_len, amp_b_len = jaccard_sim_seq(st_a, st_b)
        jaccard_bp, num_shared_bps = jaccard_sim_bp(bplist_a, bplist_b, d)
        featList.extend([jaccard_seq, jaccard_bp, num_shared_bps, len(bplist_a), len(bplist_b), amp_olap_len,
                         amp_a_len, amp_b_len])
        pcent = score_to_percentile(s, background_scores)
        pval = score_to_pval(s)
        if s != 0:
            outdata.append([a, b, s, pcent, pval] + featList)
        # outfile.write("\t".join([a, b, str(s), str(pcent), str(pval)] + [str(x) for x in featList]) + "\n")


# check if aa data repo set, construct LC datatabase
def set_lcd(AA_DATA_REPO, no_LC_filter):
    fDict = {}
    with open(AA_DATA_REPO + "file_list.txt") as infile:
        for line in infile:
            fields = line.strip().rsplit()
            fDict[fields[0]] = fields[1]

    lcPath = AA_DATA_REPO + fDict["mapability_exclude_filename"]
    cg5Path = AA_DATA_REPO + fDict["conserved_regions_filename"]

    lcD = defaultdict(IntervalTree)
    cg5D = defaultdict(IntervalTree)
    if not no_LC_filter:
        lcD = buildLCDatabase(lcPath)
        cg5D = build_CG5_database(cg5Path)

    return lcD, cg5D

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute similarity between overlapping AA amplicons")
    parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, or GRCh38",
                        choices=["hg19", "GRCh37", "hg38", "GRCh38", "mm10", "GRCm38"], required=True)
    # parser.add_argument("--min_size", type=float, help="Minimum cycle size (in bp) to consider as valid amplicon",
    #                     default=5000)
    parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: \
    samplename /path/to/sample_amplicon1_cycles.txt /path/to/sample_amplicon1_graph.txt", required=True)
    parser.add_argument("-o", help="Output filename prefix")
    parser.add_argument("--min_cn", type=float, help="Minimum CN to consider as amplification (default 4.5",
                        default=4.5)
    # parser.add_argument("-p","--pairs", help="Amplicons to compute similarity for. Each line formatted as: \
    # samplename_amplicon1 samplename_amplicon2 ...", required=True)
    parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files",
                        action='store_true', default=False)
    parser.add_argument("--no_LC_filter", help="Do not filter low-complexity cycles. Not recommended to set this flag.",
                        action='store_true', default=False)
    parser.add_argument("--min_de", type=int, help="Set the minimum number of discordant edges in the amplicon "
                                                   "required to be considered for similarity (default 1).", default=1)
    parser.add_argument("--subset_bed", help="Restrict the similarity calculation to the regions in the bed file "
                                             "provided to this argument.")
    # parser.add_argument("--cycle_similarity", help="Similarity calculations for the similarity of paths/cycles in AA "
    #                                                "in cycles files", action='store_true', default=False)
    parser.add_argument("--include_path_in_amplicon_name", help="Include path of file when reporting amplicon name",
                        action='store_true', default=False)
    parser.add_argument("--classification_file", help="Path to amplicon_classification_profiles.tsv file")
    parser.add_argument("--required_classifications", help="Amplicons must have received one or more or the following "
                        "classifications. Requires --classification_file.", choices=["ecDNA", "BFB",
                        "Complex non-cyclic", "Linear amplification", "any"], nargs='+')
    parser.add_argument("--required_amplicon_ID", help="Require that any similarity scores reported are including this "
                                                       "amplicon ID", default=None)

    args = parser.parse_args()
    if args.required_classifications and not args.classification_file:
        sys.stderr.write("--classification_file must be present with --require_classifications\n")
        sys.exit(1)
    elif args.required_classifications:
        if args.required_classifications == ["any", ]:
            required_classes = {"ecDNA", "BFB", "Complex non-cyclic", "Linear amplification"}
        else:
            required_classes = set(args.required_classifications)

    else:
        required_classes = set()

    if args.classification_file:
        a2class = read_classifications(args.classification_file)
    else:
        a2class = defaultdict(set)

    add_chr_tag = args.add_chr_tag
    cn_cut = args.min_cn
    min_de = args.min_de
    region_subset_ivald = {}

    if not args.o:
        args.o = os.path.basename(args.input).rsplit(".")[0]

    if args.ref == "hg38": args.ref = "GRCh38"

    # check if aa data repo set, construct LC datatabase
    try:
        AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + args.ref + "/"

    except KeyError:
        sys.stderr.write("$AA_DATA_REPO not set. Please see AA installation instructions.\n")
        sys.exit(1)

    lcD, cg5D = set_lcd(AA_DATA_REPO, args.no_LC_filter)

    if args.subset_bed:
        region_subset_ivald = bed_to_interval_dict(args.subset_bed, add_chr_tag)

    s2a_graph = readFlist(args.input, region_subset_ivald, lcD, cg5D, cn_cut, args.include_path_in_amplicon_name,
                          required_classes, a2class)

    pairs = get_pairs(s2a_graph, args.required_amplicon_ID)

    with open(args.o + "_amplicon_similarity_scores.tsv", 'w') as outfile:
        # [as1, as2, nS1, nS2, bpd1, bpd2]
        outfile.write("Amp1\tAmp2\tSimilarityScore\tSimScorePercentile\tSimScorePvalue\tAsymmetricScore1\t"
                      "AsymmetricScore2\tGenomicSegmentScore1\tGenomicSegmentScore2\tBreakpointScore1\t"
                      "BreakpointScore2\tJaccardGenomicSegment\tJaccardBreakpoint\tNumSharedBPs\tAmp1NumBPs\t"
                      "Amp2NumBPs\tAmpOverlapLen\tAmp1AmpLen\tAmp2AmpLen\n")
        outdata = []
        compute_similarity(outfile, s2a_graph, pairs, outdata)
        outdata.sort(key=lambda x: (x[2], x[1], x[0]), reverse=True)
        for l in outdata:
            outfile.write("\t".join([str(x) for x in l]) + "\n")
