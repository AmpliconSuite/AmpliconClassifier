from collections import defaultdict
import os
import warnings

from intervaltree import IntervalTree, Interval

lookup = str.maketrans("ACGTRYKM", "TGCAYRMK")


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

        return "{}:{}{} | {}:{}{}\t{}\t{}".format(
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

                if cn > 7:
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


# determine if a cyclic path front and back are connected by everted edge
# return true if connected, false if not connected or not a cyclic path
def amp_encompassed(cycle, segSeqD, graphf, add_chr_tag):
    # min seg and max seg in cycle -
    absSegs = set(abs(x) for x in cycle)
    minSeg, maxSeg = min(absSegs), max(absSegs)
    if minSeg == 0:
        return False

    # are they connected in same orientation?
    ca, pal = segSeqD[minSeg][0], segSeqD[minSeg][1]
    cb, pbr = segSeqD[maxSeg][0], segSeqD[maxSeg][2]

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
                    print("bridged a gap in cycle using known database: " + str(repCyc))
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
                        print("bridged a gap in cycle: " + str(repCyc))
                        return repCyc

                elif i > 0 and j > 0:
                    chri, posi = segSeqD[abs(i)][0], segSeqD[abs(i)][2]
                    chrj, posj = segSeqD[abs(j)][0], segSeqD[abs(j)][1]
                    if xt[0] == chrj and abs(posj - xt[1]) < 1000 and yt[0] == chri and abs(posi - yt[1]):
                        repCyc = repCyc[1:-1]
                        print("bridged a gap in cycle: " + str(repCyc))
                        return repCyc

    return repCyc


def parseCycle(cyclef, graphf, add_chr_tag, lcD, patch_links):
    xt, yt, ee_spans = get_amp_outside_bounds(graphf, add_chr_tag)
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

                currCycle = repair_cycle(num_ss, segSeqD, patch_links, xt, yt, ee_spans)
                uid = ss + "," + cd["Copy_count"]
                if uid in seenCycs:
                    print(cyclef + " duplicate cycle encountered")

                else:
                    cycleList.append(currCycle)
                    seenCycs.add(uid)
                    cycleCNs.append(float(cd["Copy_count"]))

    return segSeqD, cycleList, cycleCNs


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
                    print("Bad formatting in: ", line)
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

