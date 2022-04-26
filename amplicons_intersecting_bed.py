#!/usr/bin/env python3

import argparse
from collections import defaultdict

from intervaltree import IntervalTree


def readFlist(filelist):
    flist = []
    with open(filelist) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit()
                if len(fields) < 2 or len(fields) > 3:
                    print("Bad formatting in: ", line)
                else:
                    flist.append(fields)

    return flist


def intervals_from_graph(gfile):
    ivaltd = defaultdict(IntervalTree)
    with open(gfile) as infile:
        for line in infile:
            if line.startswith('sequence'):
                fields = line.rsplit()
                c1, p1 = fields[1].rsplit(':')
                p1 = int(p1[:-1])
                c2, p2 = fields[2].rsplit(':')
                p2 = int(p2[:-1]) + 1

                if not c1.startswith('chr'):
                    c1 = 'chr' + c1

                ivaltd[c1].addi(p1, p2)

    return ivaltd


def bed_to_ivals(bedf):
    ivaltd = defaultdict(IntervalTree)
    with open(bedf) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            if not fields[0].startswith('chr'):
                fields[0] = 'chr' + fields[0]

            ivaltd[fields[0]].addi(int(fields[1]), int(fields[2]))

    return ivaltd


def has_overlap(ivaltd1, ivaltd2):
    for k, ivalt1 in ivaltd1.items():
        ivalt2 = ivaltd2[k]
        for ival in ivalt1:
            if ivalt2[ival.begin:ival.end]:
                return True

    return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subsect all amplicons overlapping one or more regions from bed")
    parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: "
                        "sample_name cycles.txt graph.txt", required=True)
    parser.add_argument("-b", "--bed", help="Bed file of regions amplicon must have at least one overlap with",
                        required=True)

    args = parser.parse_args()

    bed_ivalt = bed_to_ivals(args.bed)
    flist = readFlist(args.input)
    for fs in flist:
        graph_ivalt = intervals_from_graph(fs[2])
        if has_overlap(bed_ivalt, graph_ivalt):
            print("\t".join(fs))
