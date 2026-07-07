from collections import defaultdict
import logging

logger = logging.getLogger(__name__)
import os
import warnings

from intervaltree import IntervalTree, Interval

# Module-level globals
lookup = str.maketrans("ACGTRYKM", "TGCAYRMK")


class ConfigVars:
        tot_min_del = 5000
        minCycleSize = 5000
        compCycContCut = 50000
        anyCycContcut = 10000
        ampLenOverMinCN = 5000
        cycCut = 0.12
        compCut = 0.3
        min_amp_cn = 4.5
        sig_amp = 7
        high_amp = 12
        max_segdup_size = 1000000
        segdup_max_extra_fraction = 0.25
        decomposition_strictness = 0.1
        lc_cycle_max_bp_fraction = 0.10
        lc_cycle_max_breakend_fraction = 0.10
        # bfb-related items
        min_fb_read_prop = 0.25
        fb_break_weight_prop = 0.3
        fb_dist_cut = 50000
        min_flow = 1.0
        max_foldback_edges_qc = 100
        min_bp_per_foldback = 50000
        max_nonbfb_break_weight = 0.5
        min_bfb_cycle_weight_ratio = 0.6
        bfbarchitect_max_score = 2.8
        bfbarchitect_cycle_overlap_threshold = 0.95
        bfbarchitect_min_lp_bound = 25
        bfbarchitect_max_region_segments = 50
        bfbarchitect_whole_graph_max_sequence_edges = 150
        bfbarchitect_whole_graph_max_discordant_edges = 50
        bfbarchitect_whole_graph_min_foldback_fraction = 0.1
        bfbarchitect_whole_graph_max_large_chrom_count = 3
        bfbarchitect_whole_graph_min_large_chrom_bp = 250000

        # FAN (Focal amplification in neochromosome) related items
        fan_ecDNA_min_cn = 60

        # TID-related items
        tid_cn_ratio_max = 2.5
        tid_cn_ratio_max_tight_fb = 4.0
        tid_tight_fb_size = 50000
        tid_max_foldback_span = 1000000
        tid_close_endpoint_dist = 25000
        tid_max_inner_cn = 12


class breakpoint(object):
    def __init__(self, lchrom, lpos, rchrom, rpos, cn=None, support=None, ldir=None, rdir=None, homlen="NA", homseq="NA"):
        # Use __setattr__ before setting _hash to avoid warnings during initialization
        object.__setattr__(self, '_hash', None)
        self.lchrom = lchrom
        self.lpos = lpos
        self.rchrom = rchrom
        self.rpos = rpos
        self.cn = cn
        self.support = support
        self.ldir = ldir
        self.rdir = rdir
        self.homlen = homlen
        self.homseq = homseq

    def __setattr__(self, name, value):
        if name != '_hash' and hasattr(self, '_hash') and self._hash is not None:
            warnings.warn(
                "Modifying attribute '{}' after object has been hashed. "
                "This may cause unexpected behavior if the object is being used in "
                "sets or as a dictionary key.".format(name),
                RuntimeWarning,
                stacklevel=2
            )
            # Clear the cached hash
            object.__setattr__(self, '_hash', None)
        object.__setattr__(self, name, value)

    def __eq__(self, other):
        if not isinstance(other, breakpoint):
            return NotImplemented
        return (
            self.lchrom == other.lchrom and
            self.lpos == other.lpos and
            self.rchrom == other.rchrom and
            self.rpos == other.rpos and
            self.cn == other.cn and
            self.support == other.support and
            self.ldir == other.ldir and
            self.rdir == other.rdir and
            self.homlen == other.homlen and
            self.homseq == other.homseq
        )

    def __hash__(self):
        if self._hash is None:
            self._hash = hash((
                self.lchrom,
                self.lpos,
                self.rchrom,
                self.rpos,
                self.cn,
                self.support,
                self.ldir,
                self.rdir,
                self.homlen,
                self.homseq
            ))
        return self._hash

    def __str__(self):
        ldir = self.ldir if self.ldir is not None else ""
        rdir = self.rdir if self.rdir is not None else ""
        cn = str(self.cn) if self.cn is not None else "NA"
        support = str(self.support) if self.support is not None else "NA"

        return "{}:{}{} | {}:{}{} | {} | {}".format(
            self.lchrom, self.lpos, ldir, self.rchrom, self.rpos, rdir, cn, support)

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

# ------------------------------------------------------------


class amped_gene(object):
    def __init__(self, name, ncbi_id, has5p_end, has3p_end, gene_cn, is_oncogene):
        self.name = name
        self.ncbi_id = ncbi_id
        self.has5p_end = has5p_end
        self.has3p_end = has3p_end
        self.gene_cn = gene_cn
        self.is_oncogene = is_oncogene
        self.trunc_set = set()

    def get_truncation_status(self):
        tset = []
        if not self.has5p_end:
            tset.append("5p")
        if not self.has3p_end:
            tset.append("3p")

        if tset:
            tstring = "_".join(tset)
        else:
            tstring = "None"

        return tstring


def merge_intervals(feature_dict, tol=1):
    for item, usort_intd in feature_dict.items():
        for chrom, usort_ints in usort_intd.items():
            if not usort_ints:
                continue

            # sort ints
            sort_ints = sorted(usort_ints)
            # merge sorted ints
            mi = [sort_ints[0]]
            for ival in sort_ints[1:]:
                if ival[0] <= mi[-1][1] + tol:
                    ui = (mi[-1][0], max(ival[1], mi[-1][1]))
                    mi[-1] = ui

                else:
                    mi.append(ival)

            feature_dict[item][chrom] = mi


def is_human_viral_hybrid(ref, cycle, segSeqD):
    chromList = [segSeqD[abs(ind)][0] for ind in cycle if ind != 0]
    return ref == "GRCh38_viral" and any([not x.startswith("chr") for x in chromList]) and any([x.startswith("chr") for x in chromList])


def is_viral(ref, cycle, segSeqD):
    chromList = [segSeqD[abs(ind)][0] for ind in cycle if ind != 0]
    return ref == "GRCh38_viral" and any([not x.startswith("chr") for x in chromList])


def isCircular(cycle):
    return cycle[0] != 0 or cycle[-1] != 0


def get_amp_outside_bounds(graphf, add_chr_tag):
    # get the interval of the first amp (x)
    # get the interval of the last amp (y)
    # is there an everted edge linking x to y?
    xc, xs = None, 0
    yc, ye = None, 0
    ee_spans = False
    with open(graphf) as infile:
        for line in infile:
            fields = line.rsplit()
            if line.startswith("sequence"):
                c, s, e = fields[1].rsplit(":")[0], int(fields[1].rsplit(":")[1][:-1]), int(fields[2].rsplit(":")[1][:-1])+1
                cn = float(fields[3])
                if add_chr_tag and not c.startswith('chr'):
                    c = "chr" + c

                if cn > ConfigVars.sig_amp:
                    if not xc:
                        xc, xs = c, s

                    yc, ye = c, e

            elif line.startswith("discordant") and not xc is None and not yc is None:
                s1, s2 = fields[1].rsplit("->")
                c1, p1, d1 = s1.rsplit(":")[0],  int(s1.rsplit(":")[1][:-1]), s1.rsplit(":")[1][-1]
                c2, p2, d2 = s2.rsplit(":")[0], int(s2.rsplit(":")[1][:-1]), s2.rsplit(":")[1][-1]
                if add_chr_tag and not c1.startswith('chr'):
                    c1 = "chr" + c1
                    c2 = "chr" + c2

                # first must be +, second must be - for everted orientation
                if d1 == "+" and d2 == "-" and c1 == c2 and p1 > p2 and c1 == yc and c2 == xc:
                    # first must fall within 1kbp of yc,ye
                    # second must fall within 1kbp of xc,xs
                    if ye - 1000 < p1 < ye + 1000 and xs - 1000 < p2 < xs + 1000:
                        ee_spans = True

    return (xc, xs), (yc, ye), ee_spans


def _parse_graph_vertex(vertex, add_chr_tag=False):
    chrom, pos_dir = vertex.rsplit(":", 1)
    if add_chr_tag and not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return chrom, int(pos_dir[:-1]), pos_dir[-1]


def _discordant_edge_key(v1, v2):
    return tuple(sorted((v1, v2)))


def get_discordant_edge_keys(graphf, add_chr_tag=False):
    discordant_edges = set()
    with open(graphf) as infile:
        for line in infile:
            if not line.startswith("discordant"):
                continue
            fields = line.rstrip().rsplit()
            if len(fields) < 2:
                continue
            try:
                v1, v2 = fields[1].rsplit("->", 1)
                discordant_edges.add(_discordant_edge_key(
                    _parse_graph_vertex(v1, add_chr_tag),
                    _parse_graph_vertex(v2, add_chr_tag),
                ))
            except (IndexError, ValueError):
                continue

    return discordant_edges


def _signed_segment_entry_vertex(signed_seg, segSeqD):
    chrom, start, end = segSeqD[abs(signed_seg)]
    if signed_seg > 0:
        return chrom, start, "-"
    return chrom, end, "+"


def _signed_segment_exit_vertex(signed_seg, segSeqD):
    chrom, start, end = segSeqD[abs(signed_seg)]
    if signed_seg > 0:
        return chrom, end, "+"
    return chrom, start, "-"


def _walk_junction_pairs(num_ss):
    pairs = list(zip(num_ss[:-1], num_ss[1:]))
    if 0 not in num_ss and len(num_ss) > 1:
        pairs.append((num_ss[-1], num_ss[0]))
    return pairs


def _lc_overlap_bp(lcD, chrom, start, end):
    return sum(max(0, min(end, iv.end) - max(start, iv.begin)) for iv in lcD[chrom].overlap(start, end))


def lc_cycle_content(num_ss, segSeqD, lcD, discordant_edge_keys):
    total_bp = 0
    lc_bp = 0
    first_lc_overlap = None
    for signed_seg in num_ss:
        if signed_seg == 0:
            continue
        chrom, start, end = segSeqD[abs(signed_seg)]
        total_bp += max(0, end - start)
        overlap_bp = _lc_overlap_bp(lcD, chrom, start, end)
        if overlap_bp:
            lc_bp += overlap_bp
            if first_lc_overlap is None:
                first_lc_overlap = (chrom, start, end)

    total_breakends = 0
    lc_breakends = 0
    for a, b in _walk_junction_pairs(num_ss):
        if a == 0 or b == 0:
            continue
        v1 = _signed_segment_exit_vertex(a, segSeqD)
        v2 = _signed_segment_entry_vertex(b, segSeqD)
        if _discordant_edge_key(v1, v2) not in discordant_edge_keys:
            continue

        for chrom, pos, _strand in (v1, v2):
            total_breakends += 1
            if lcD[chrom].overlaps(pos, pos + 1):
                lc_breakends += 1

    bp_fraction = lc_bp / float(total_bp) if total_bp else 0.0
    breakend_fraction = lc_breakends / float(total_breakends) if total_breakends else 0.0
    return {
        "total_bp": total_bp,
        "lc_bp": lc_bp,
        "bp_fraction": bp_fraction,
        "total_breakends": total_breakends,
        "lc_breakends": lc_breakends,
        "breakend_fraction": breakend_fraction,
        "first_lc_overlap": first_lc_overlap,
    }


def get_diff(e1, e2, segSeqD):
    p1_abs = segSeqD[abs(e1)]
    p2_abs = segSeqD[abs(e2)]
    if e1 == 0 or e2 == 0:
        return 1

    if p1_abs[0] != p2_abs[0]:
        return None

    p1_end = p1_abs[2] if e1 > 0 else p1_abs[1]
    p2_start = p2_abs[1] if e2 > 0 else p2_abs[2]
    return abs(p2_start - p1_end)


def pos_lies_on_cn_change(chrom, pos, graph_cns, min_delta=2, tol=150):
    hits = graph_cns[chrom][pos-tol:pos+tol+1]
    if not hits:
        warnings.warn("Attempted to check boundary on CN change, but region not found in graph.")
        return False

    cns = [hit.data for hit in hits if hit.end - hit.begin > tol]
    return max(cns) - min(cns) >= min_delta


def min_max_cycle_posns(cycle, segSeqD):
    # min seg and max seg in cycle -
    absSegs = set(abs(x) for x in cycle if x != 0)
    minSeg, maxSeg = min(absSegs), max(absSegs)

    ca, pal = segSeqD[minSeg][0], segSeqD[minSeg][1]
    cb, pbr = segSeqD[maxSeg][0], segSeqD[maxSeg][2]

    return minSeg, maxSeg, (ca, pal), (cb, pbr)


# determine if a cyclic path front and back are connected by everted edge
# return true if connected, false if not connected or not a cyclic path
def amp_encompassed(cycle, segSeqD, graphf, add_chr_tag):
    if not isCircular(cycle):
        return False

    minSeg, maxSeg, (ca, pal), (cb, pbr) = min_max_cycle_posns(cycle, segSeqD)
    # is there an everted edge joining them?
    with open(graphf) as infile:
        for line in infile:
            fields = line.rsplit()
            if line.startswith("discordant"):
                s1, s2 = fields[1].rsplit("->")
                c1, p1, d1 = s1.rsplit(":")[0],  int(s1.rsplit(":")[1][:-1]), s1.rsplit(":")[1][-1]
                c2, p2, d2 = s2.rsplit(":")[0], int(s2.rsplit(":")[1][:-1]), s2.rsplit(":")[1][-1]
                if add_chr_tag and not c1.startswith('chr'):
                    c1 = "chr" + c1
                    c2 = "chr" + c2

                # first must be +, second must be -
                if d1 == "+" and d2 == "-" and ca == c2 and cb == c1:
                    if pbr - 1 < p1 < pbr + 1 and pal - 1 < p2 < pal + 1:
                        return True

    return False


# Input and data prep
def bpgEdgeToCycles(bp, posCycleLookup):
    lCycles = set([x.data for x in posCycleLookup[bp.lchrom][bp.lpos]])
    rCycles = set([x.data for x in posCycleLookup[bp.rchrom][bp.rpos]])
    return lCycles, rCycles


def rotate_cycle_to_largest_jump(cycle, segSeqD):
    if not isCircular(cycle) or len(cycle) == 1:
        return cycle

    diffs = [get_diff(cycle[-1], cycle[0], segSeqD)]
    for ind, seg in enumerate(cycle[:-1]):
        a_b_diff = get_diff(cycle[ind], cycle[ind + 1], segSeqD)
        diffs.append(a_b_diff)

    # Find best rotation index based on diffs
    best_index = 0
    best_score = float('-inf')

    for i, diff in enumerate(diffs):
        # Prioritize None over numbers, and then pick the largest number
        if diff is None:
            best_index = i
            break
        elif diff > best_score:
            best_score = diff
            best_index = i

    # Rotate so that cycle[best_index] becomes the first element
    rotated_cycle = cycle[best_index:] + cycle[:best_index]
    return rotated_cycle


# address formatting issues and known link misses
def repair_cycle(cycle, segSeqD, patch_links, xt, yt, ee_spans):
    repCyc = cycle

    # repair issue with interior source edges
    if 0 in cycle and cycle[0] != 0:
        zero_ind = cycle.index(0)
        repCyc = cycle[zero_ind + 1:] + cycle[:zero_ind + 1]

    # repair issues with unlinked deletion
    if len(repCyc) > 3 and repCyc[0] == 0:
        # repair small unlinked deletions from known database (patch_links)
        directional_a, directional_b = repCyc[1], repCyc[-2]
        a, b = abs(directional_a), abs(directional_b)
        if directional_a > 0:
            chra, posa = segSeqD[a][0], segSeqD[a][1]
        else:
            chra, posa = segSeqD[a][0], segSeqD[a][2]

        if directional_b > 0:
            chrb, posb = segSeqD[b][0], segSeqD[a][2]
        else:
            chrb, posb = segSeqD[b][0], segSeqD[b][1]

        spair = sorted([(chra, posa), (chrb, posb)])
        for x in patch_links:
            if x[0] == spair[0][0] and x[2] == spair[1][0]:
                if x[1].overlaps(spair[0][1]) and x[3].overlaps(spair[1][1]):
                    repCyc = repCyc[1:-1]
                    logger.info("bridged a gap in cycle using known database: " + str(repCyc))
                    return repCyc

        # now check if amplicon entirely enclosed in everted edge with high CN (suggests missing interior edge).
        if ee_spans:
            for i, j in zip(repCyc[1:-2], repCyc[2:-1]):
                # check if same direction
                if i < 0 and j < 0:
                    chri, posi = segSeqD[abs(i)][0], segSeqD[abs(i)][1]
                    chrj, posj = segSeqD[abs(j)][0], segSeqD[abs(j)][2]
                    if xt[0] == chri and abs(posi - xt[1]) < 1000 and yt[0] == chrj and abs(posj - yt[1]):
                        repCyc = repCyc[1:-1]
                        logger.info("bridged a gap in cycle: " + str(repCyc))
                        return repCyc

                elif i > 0 and j > 0:
                    chri, posi = segSeqD[abs(i)][0], segSeqD[abs(i)][2]
                    chrj, posj = segSeqD[abs(j)][0], segSeqD[abs(j)][1]
                    if xt[0] == chrj and abs(posj - xt[1]) < 1000 and yt[0] == chri and abs(posi - yt[1]):
                        repCyc = repCyc[1:-1]
                        logger.info("bridged a gap in cycle: " + str(repCyc))
                        return repCyc

    return repCyc


def parseCycle(cyclef, graphf, add_chr_tag, lcD, patch_links):
    xt, yt, ee_spans = get_amp_outside_bounds(graphf, add_chr_tag)
    discordant_edge_keys = get_discordant_edge_keys(graphf, add_chr_tag)
    segSeqD = {0: (None, 0, 0)}
    cycleList = []
    cycleCNs = []
    lc_filtered_cycles = 0
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

            elif line.startswith("Cycle=") or line.startswith("Path="):
                # CoRAL labels non-cyclic walks "Path=" and cyclic walks "Cycle="; AA labels
                # both "Cycle=". Accept either - cyclic vs. path is inferred from the 0 segment.
                # Require the '=' so CoRAL's "Path constraint" lines are not mistaken for walks.
                cf = [tuple(x.rsplit("=")) for x in line.rstrip().rsplit(";")]
                cd = dict(cf)
                ss = cd["Segments"]
                num_ss = [int(x[-1] + x[:-1]) for x in ss.rsplit(",")]
                pop_inds = []
                for seg_ind, seg in enumerate(num_ss):
                    t = segSeqD[abs(seg)]
                    if lcD[t[0]].overlaps(t[1], t[2]):
                        if num_ss[0] == 0 and (seg_ind == 1 or seg_ind == len(num_ss) - 2) and len(num_ss) > 3:
                            pop_inds.append(seg_ind)
                            continue

                pop_ind_set = set(pop_inds)
                eval_num_ss = [seg for seg_ind, seg in enumerate(num_ss) if seg_ind not in pop_ind_set]
                lc_stats = lc_cycle_content(eval_num_ss, segSeqD, lcD, discordant_edge_keys)
                lc_bp_filter = (
                    lc_stats["lc_bp"] > 0 and
                    lc_stats["bp_fraction"] >= ConfigVars.lc_cycle_max_bp_fraction
                )
                lc_breakend_filter = (
                    lc_stats["lc_breakends"] > 0 and
                    lc_stats["breakend_fraction"] >= ConfigVars.lc_cycle_max_breakend_fraction
                )

                if lc_bp_filter or lc_breakend_filter:
                    first_overlap = lc_stats["first_lc_overlap"]
                    if first_overlap:
                        overlap_msg = " | overlaps {}:{}-{}".format(
                            str(first_overlap[0]), str(first_overlap[1]), str(first_overlap[2])
                        )
                    else:
                        overlap_msg = ""
                    logger.info(
                        "Cycle was LC: {}{} | LC_bp_fraction={:.4f}; LC_breakend_fraction={:.4f}".format(
                            line.rstrip(), overlap_msg, lc_stats["bp_fraction"], lc_stats["breakend_fraction"]
                        )
                    )
                    lc_filtered_cycles += 1
                    continue

                elif pop_inds:
                    for seg_ind in pop_inds[::-1]:
                        num_ss.pop(seg_ind)

                initCycle = repair_cycle(num_ss, segSeqD, patch_links, xt, yt, ee_spans)
                currCycle = rotate_cycle_to_largest_jump(initCycle, segSeqD)
                uid = ss + "," + cd["Copy_count"]
                if uid in seenCycs:
                    logger.info(cyclef + " duplicate cycle encountered")

                else:
                    cycleList.append(currCycle)
                    seenCycs.add(uid)
                    cycleCNs.append(float(cd["Copy_count"]))

    return segSeqD, cycleList, cycleCNs, lc_filtered_cycles


def get_graph_cns(gfile, add_chr_tag):
    cns_dict = defaultdict(IntervalTree)
    with open(gfile) as infile:
        for line in infile:
            if line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                cn = float(fields[3])
                c, p1 = fields[1].rsplit(":")
                p1 = int(p1[:-1])
                c, p2 = fields[2].rsplit(":")
                p2 = int(p2[:-1]) + 1
                if add_chr_tag and not c.startswith('chr'):
                    c = 'chr' + c

                cns_dict[c].addi(p1, p2, cn)

    return cns_dict


# build a lookup for position to list of cycles hitting it
def buildPosCycleLookup(cycles, segSeqD):
    posCycleLookup = defaultdict(IntervalTree)
    for ind, c in enumerate(cycles):
        for seg in c:
            if seg != 0:
                chrom, s, e = segSeqD[abs(seg)]
                posCycleLookup[chrom][s:e + 1] = ind

    return posCycleLookup


def buildLCDatabase(mappabilityFile, filtsize=7500):
    lcD = defaultdict(IntervalTree)
    with open(mappabilityFile) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if fields:
                chrom, s, e = fields[0], int(fields[1]), int(fields[2])
                if e - s > filtsize:
                    lcD[chrom].addi(s, e)

    return lcD


def build_CG5_database(cg5_file):
    cg5D = defaultdict(IntervalTree)
    with open(cg5_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if not fields:
                continue
            chrom, s, e = fields[0], int(fields[1]), int(fields[2])
            cg5D[chrom].addi(s, e)

    return cg5D


# check if aa data repo set, construct LC datatabase
def set_lcd(AA_DATA_REPO, no_LC_filter=False):
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


def readFlist(filelist):
    flist = []
    with open(filelist) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit("\t")
                if len(fields) < 2:
                    logger.warning("Bad formatting in: ", line)
                else:
                    flist.append(fields)

    return flist


def read_patch_regions(ref):
    dp = os.path.dirname(os.path.abspath(__file__)) + "/resources/"
    patch_links = []
    with open(dp + "patch_regions.tsv") as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            if fields[0] == ref:
                patch_links.append([fields[1], Interval(int(fields[2]), int(fields[3])),
                                    fields[4], Interval(int(fields[5]), int(fields[6]))])

    return patch_links


def write_patched_graph(gfile, patch_links, outpath):
    """
    Write a copy of gfile to outpath with patch-region sequence edge CNs replaced by the
    mean of their genomic neighbors. Flanking concordant edge CNs are updated to match.
    Returns True if any patches were applied, False otherwise.
    """
    if not patch_links:
        return False

    with open(gfile) as f:
        lines = f.readlines()

    # Pass 1: collect sequence edges and identify patch segments
    segs = []  # (chrom, s, e, cn, line_index)
    in_seq = False
    for i, line in enumerate(lines):
        t = line.rstrip()
        p = t.split('\t')
        if t.startswith('SequenceEdge:'):   in_seq = True;  continue
        if t.startswith('BreakpointEdge:'): in_seq = False; continue
        if in_seq and len(p) >= 4 and p[0] == 'sequence':
            try:
                c = p[1].split(':')[0]
                pl = int(p[1].split(':')[1][:-1]); pr = int(p[2].split(':')[1][:-1])
                segs.append((c, min(pl, pr), max(pl, pr), float(p[3]), i))
            except:
                pass

    by_chrom = defaultdict(list)
    for rec in segs:
        by_chrom[rec[0]].append(rec)
    for c in by_chrom:
        by_chrom[c].sort(key=lambda x: x[1])

    cn_replacements = {}    # line_index -> imputed cn
    patch_boundary_pos = {} # (chrom, pos) -> imputed cn, for concordant edge correction

    for chrom, csegs in by_chrom.items():
        for i in range(1, len(csegs) - 1):
            c, s, e, cn, li = csegs[i]
            for x in patch_links:
                if x[0] == chrom and x[2] == chrom and x[1].overlaps(s - 1) and x[3].overlaps(e + 1):
                    imputed = (csegs[i - 1][3] + csegs[i + 1][3]) / 2.0
                    cn_replacements[li] = imputed
                    for pos in (s - 1, s, e, e + 1):
                        patch_boundary_pos[(chrom, pos)] = imputed
                    logger.info("patching graph {}:{}-{} CN {:.2f} -> {:.2f}".format(
                        chrom, s, e, cn, imputed))
                    break

    if not cn_replacements:
        return False

    # Pass 2: write corrected graph
    in_seq = in_bp = False
    with open(outpath, 'w') as out:
        for i, line in enumerate(lines):
            t = line.rstrip()
            p = t.split('\t')
            if t.startswith('SequenceEdge:'):   in_seq = True;  in_bp = False
            elif t.startswith('BreakpointEdge:'): in_seq = False; in_bp = True

            if in_seq and i in cn_replacements:
                p[3] = repr(cn_replacements[i])
                out.write('\t'.join(p) + '\n')
            elif in_bp and p[0] == 'concordant' and len(p) >= 3:
                try:
                    v1, v2 = p[1].split('->')
                    c1 = v1.split(':')[0]; pos1 = int(v1.split(':')[1][:-1])
                    new_cn = patch_boundary_pos.get((c1, pos1))
                    if new_cn is None:
                        c2 = v2.split(':')[0]; pos2 = int(v2.split(':')[1][:-1])
                        new_cn = patch_boundary_pos.get((c2, pos2))
                    if new_cn is not None:
                        p[2] = repr(new_cn)
                        out.write('\t'.join(p) + '\n')
                        continue
                except:
                    pass
                out.write(line)
            else:
                out.write(line)

    return True


def get_ncrna_file_loc(ref):
    dp = os.path.dirname(os.path.abspath(__file__)) + "/resources/"
    simple_rname = ref
    if ref == "GRCh37":
        simple_rname = 'hg19'
    if ref == "GRCh38_viral":
        simple_rname = 'GRCh38'

    fname = "{}gencode_{}_long_noncoding_RNAs.gff3.gz".format(dp, simple_rname)
    return fname


def set_config_vars(config, args):
    for key in config:
        if hasattr(ConfigVars, key):
            setattr(ConfigVars, key, config[key])

    # Optional overrides from command line args
    if hasattr(args, 'decomposition_strictness') and args.decomposition_strictness is not None:
        ConfigVars.decomposition_strictness = args.decomposition_strictness

    if hasattr(args, 'min_size') and args.min_size is not None:
        ConfigVars.minCycleSize = args.min_size

    if hasattr(args, 'min_flow') and args.min_flow is not None:
        ConfigVars.min_flow = args.min_flow
