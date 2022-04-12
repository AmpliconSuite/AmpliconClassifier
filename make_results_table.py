#!/usr/bin/env python3

import argparse
import bisect
from collections import defaultdict
import os
import sys

from intervaltree import IntervalTree


# The input file must be the source of the classification file calls.

parser = argparse.ArgumentParser(description="Organize AC results into a table")
parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: "
                    "sample_name cycles.txt graph.txt.", required=True)
parser.add_argument("--classification_file", help="Path to amplicon_classification_profiles.tsv file", required=True)
args = parser.parse_args()

output_head = ["Sample name", "AA amplicon number", "Feature ID", "Classification", "Location", "AA PNG file",
               "AA PDF file"]

output_table_lines = [output_head, ]
with open(args.input) as input_file, open(args.classification_file) as classification_file:
    classBase = args.classification_file.rsplit("_amplicon_classification_profiles.tsv")[0]
    classBedDir = classBase + "_classification_bed_files/"
    class_head = next(classification_file).rstrip().rsplit("\t")
    for input_line, classification_line in zip(input_file, classification_file):
        input_fields = input_line.rstrip().rsplit()
        sample_name = input_fields[0].rsplit("_amplicon")[0]
        amplicon_prefix = input_fields[1].rsplit("_cycles.txt")[0]
        AA_amplicon_number = amplicon_prefix.rsplit("_amplicon")[-1]

        classD = dict(zip(class_head, classification_line.rstrip().rsplit("\t")))
        ampliconID = "_".join([classD["sample_name"], classD["amplicon_number"]])
        if sample_name + "_amplicon" + AA_amplicon_number != ampliconID:
            sys.stderr.write(sample_name + "_amplicon" + AA_amplicon_number + " | " + ampliconID + "\n")
            sys.stderr.write("File ordering in " + args.input + " does not match order in "
                             + args.classification_file + "\n")
            sys.exit(1)

        # look up image locations
        # png, pdf, .... others?
        AA_png_loc = amplicon_prefix + ".png"
        AA_pdf_loc = amplicon_prefix + ".pdf"
        image_locs = [AA_png_loc, AA_pdf_loc]
        for ind, f in enumerate(image_locs):
            if not os.path.exists(f):
                sys.stderr.write("Warning: image file " + f + " not found!\n")
                image_locs[ind] = "Not found"


        amps_classes = []
        if classD["ecDNA+"] == "Positive":
            amps_classes.append(("ecDNA", int(classD["ecDNA_amplicons"])))

        if classD["BFB+"] == "Positive":
            amps_classes.append(("BFB", 1))

        elif not amps_classes and classD["amplicon_decomposition_class"] != "No amp/Invalid":
            amps_classes.append((classD["amplicon_decomposition_class"], 1))

        # Get the AC intervals
        featureData = []
        for feature, namps in amps_classes:
            for i in range(namps):
                featureID = "_".join([ampliconID, feature, str(i+1)])
                featureBed = classBedDir + featureID + "_intervals.bed"
                if not os.path.exists(featureBed):
                    sys.stderr.write("Warning: image file " + f + " not found!\n")
                    intervals = "Interval file not found"

                else:
                    interval_list = []
                    with open(featureBed) as bedfile:
                        for l in bedfile:
                            bfields = l.rstrip().rsplit("\t")
                            interval_list.append(bfields[0] + ":" + bfields[1] + "-" + bfields[2])

                    intervals = "|".join(interval_list)

                featureData.append([featureID, feature, intervals])

        for ft in featureData:
            output_table_lines.append([sample_name, AA_amplicon_number] + ft + image_locs)

with open(classBase + "_result_table.tsv", 'w') as outfile:
    for ll in output_table_lines:
        oline = "\t".join(ll) + "\n"
        outfile.write(oline)
