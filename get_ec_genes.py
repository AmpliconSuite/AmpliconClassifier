#!/usr/bin/env python3

import argparse
import sys

import pandas as pd

minCycleSize = 10000


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get genes carried on ecDNA+ samples")
    parser.add_argument("-c", "--cycles", help="AA-formatted cycles file")
    parser.add_argument("-g", "--graph", help="AA-formatted graph file")
    parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, or GRCh38",
                        choices=["hg19", "GRCh37", "hg38", "GRCh38"], required=True)
    parser.add_argument("--min_cn_flow", type=float, help="Minimum CN flow to consider as amplification", default=1)
    parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: \
    sample_name cycles.txt graph.txt")

    # parser.add_argument("--add_chr_tag", help="Add \'chr\' to the beginning of chromosome names in input files",
    #                     action='store_true')

    args = parser.parse_args()

    if not (args.cycles and args.graph) and not args.input:
        print("Need to specify (--cycles & --graph) or --input\n")
        sys.exit(1)

    if args.ref == "hg38": args.ref = "GRCh38"

    # check if aa data repo set, construct LC datatabase
    try:
        AA_DATA_REPO = os.environ["AA_DATA_REPO"] + "/" + args.ref + "/"
        fDict = {}
        with open(AA_DATA_REPO + "file_list.txt") as infile:
            for line in infile:
                fields = line.strip().rsplit()
                fDict[fields[0]] = fields[1]

        lcPath = AA_DATA_REPO + fDict["mapability_exclude_filename"]
        lcD = buildLCDatabase(lcPath)

    except KeyError:
        sys.stderr.write("$AA_DATA_REPO not set. Please see AA installation instructions.\n")
        sys.exit(1)

    if not args.input:
        tempName = args.cycles.rsplit("/")[-1].rsplit(".")[0]
        flist = [[tempName, args.cycles, args.graph]]
        if not args.o:
            args.o = os.path.basename(args.cycles).rsplit("_cycles.txt")[0]

    else:
        flist = readFlist(args.input)
        if not args.o:
            args.o = os.path.basename(args.input).rsplit(".")[0]

    minCycleSize = args.min_size