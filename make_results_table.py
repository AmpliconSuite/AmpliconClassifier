#!/usr/bin/env python3

import argparse
from collections import defaultdict
import json
import os
import shutil
import sys


def read_amplicon_gene_list(gene_file):
    amplicon_gene_dict = defaultdict(list)
    with open(gene_file) as infile:
        h = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            fd = dict(zip(h, fields))
            featureID = "_".join(fields[:3])
            if "5p" not in fd["truncated"]:
                amplicon_gene_dict[featureID].append((fd['gene'], fd['gene_cn'], eval(fd['is_canonical_oncogene'])))

    return amplicon_gene_dict


def read_complexity_scores(entropy_file):
    amplicon_complexity_dict = defaultdict(lambda: "NA")
    with open(entropy_file) as infile:
        h = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            # fd = dict(zip(h, fields))
            featureID = "_".join(fields[:3])
            amplicon_complexity_dict[featureID] = fields[4]

    return amplicon_complexity_dict


def copy_AA_files(ll):
    ldir = "files/"
    if not os.path.exists(ldir):
        os.makedirs(ldir)

    for i in range(-4, 0):
        s = ll[i]
        if s != "Not found":
            shutil.copy(s, ldir)
            ll[i] = ldir + os.path.basename(ll[i])


def html_table(output_table_lines, html_ofname):
    with open(html_ofname, 'w') as outfile:
        outfile.write("<style>\ntable, th, td {\n    border: 1px solid black;\n}\n</style>\n")
        outfile.write('<table>\n')
        for ind, ll in enumerate(output_table_lines):
            # hll = [x.replace("/opt/gpbeta_2/gp_home/", "files/") for x in ll]
            if ind != 0:
                copy_AA_files(ll)
                for i in range(-3, 0):
                    s = ll[i]
                    if s != "Not found":
                        ll[i] = "<a href=" + s + ">File</a>"
            outfile.write('<tr><td>')
            outfile.write('</td>\n    <td>'.join(ll))
            outfile.write('</td></tr>\n')

        outfile.write('</table>\n')


if __name__ == "__main__":
    # The input file must be the source of the classification file calls.
    parser = argparse.ArgumentParser(description="Organize AC results into a table")
    parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: "
                        "sample_name cycles.txt graph.txt.", required=True)
    parser.add_argument("--classification_file", help="Path to amplicon_classification_profiles.tsv file", required=True)
    parser.add_argument("--metadata_dict", help="Path to [sample]_run_metadata.json file", default="")
    args = parser.parse_args()

    output_head = ["Sample name", "AA amplicon number", "Feature ID", "Classification", "Location", "Oncogenes",
                   "Complexity score", "Reference version", "Feature BED file", "AA PNG file", "AA PDF file",
                   "Run metadata JSON"]

    output_table_lines = [output_head, ]
    with open(args.input) as input_file, open(args.classification_file) as classification_file:
        classBase = args.classification_file.rsplit("_amplicon_classification_profiles.tsv")[0]
        classBedDir = classBase + "_classification_bed_files/"
        gene_file = classBase + "_gene_list.tsv"
        entropy_file = classBase + "_feature_entropy.tsv"
        amplicon_gene_dict = read_amplicon_gene_list(gene_file)
        amplicon_complexity_dict = read_complexity_scores(entropy_file)
        if args.metadata_dict:
            metadata_dict = json.load(open(args.metadata_dict, 'r'))

        else:
            args.metadata_dict = "Not found"
            metadata_dict = defaultdict(lambda: "NA")

        class_head = next(classification_file).rstrip().rsplit("\t")
        for input_line, classification_line in zip(input_file, classification_file):
            input_fields = input_line.rstrip().rsplit()
            sample_name = input_fields[0].rsplit("_amplicon")[0]
            amplicon_prefix = input_fields[1].rsplit("_cycles.txt")[0]
            if not ":" in amplicon_prefix: amplicon_prefix.replace("//","/")
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

            # Get the AC intervals, genes and complexity
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

                    raw_glist = amplicon_gene_dict[featureID]
                    oncogenes = "|".join(sorted([g[0] for g in raw_glist if g[2]]))
                    complexity = amplicon_complexity_dict[featureID]

                    featureData.append([featureID, feature, intervals, oncogenes, complexity,
                                        metadata_dict["ref_genome"], os.path.abspath(featureBed)])

            for ft in featureData:
                output_table_lines.append([sample_name, AA_amplicon_number] + ft + image_locs + [args.metadata_dict,])

    tsv_ofname = classBase + "_result_table.tsv"
    # html_ofname = classBase + "_GenePatternNotebook_result_table.html"
    html_ofname = "index.html"

    with open(tsv_ofname, 'w') as outfile:
        for ll in output_table_lines:

            oline = "\t".join(ll) + "\n"
            outfile.write(oline)

    html_table(output_table_lines, html_ofname)
