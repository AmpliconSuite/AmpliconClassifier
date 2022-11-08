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


def read_basic_stats(basic_stats_file):
    basic_stats_dict = defaultdict(lambda: ["NA", "NA", "NA"])
    with open(basic_stats_file) as infile:
        h = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            featureID = fields[0]
            basic_stats_dict[featureID] = fields[1:]

    return basic_stats_dict


def copy_AA_files(ll, ldir):
    for i in range(-5, 0):
        s = ll[i]
        if not s.endswith("Not found") and not s.endswith("Not provided"):
            if not os.path.exists(ldir + os.path.basename(s)):
                shutil.copy(s, ldir)

            ll[i] = ldir + os.path.basename(ll[i])


def write_html_table(output_table_lines, html_ofname):
    with open(html_ofname, 'w') as outfile:
        outfile.write("<style>\ntable, th, td {\n    border: 1px solid black;\n}\n</style>\n")
        outfile.write('<table>\n')
        for ind, ll in enumerate(output_table_lines):
            # hll = [x.replace("/opt/gpbeta_2/gp_home/", "files/") for x in ll]
            if ind != 0:
                for i in range(-5, 0):
                    s = ll[i]
                    if not s.endswith("Not found") and not s.endswith("Not provided"):
                        ll[i] = "<a href=" + s + ">File</a>"
            outfile.write('<tr><td>')
            outfile.write('</td>\n    <td>'.join(ll))
            outfile.write('</td></tr>\n')

        outfile.write('</table>\n')


def write_json_dict(output_table_lines, json_ofname):
    # zip the head to the line to make a dict and dump them.
    dlist = []
    h = output_table_lines[0]
    for ll in output_table_lines[1:]:

        td = dict(zip(h, ll))
        dlist.append(td)

    with open(json_ofname, 'w') as outfile:
        json.dump(dlist, outfile, sort_keys=True)

    pass


if __name__ == "__main__":
    # The input file must be the source of the classification file calls.
    parser = argparse.ArgumentParser(description="Organize AC results into a table")
    parser.add_argument("-i", "--input", help="Path to list of files to use. Each line formatted as: "
                        "sample_name cycles.txt graph.txt.", required=True)
    parser.add_argument("--classification_file", help="Path of amplicon_classification_profiles.tsv file",
                        required=True)
    parser.add_argument("--run_metadata_file", help="Path of run metadata, [sample]_run_metadata.json file (for single"
                                                    " sample).", default="")
    parser.add_argument("--cnv_bed", help="Path of the CNV_CALLS.bed file used for this run (for single sample).",
                        default="")
    parser.add_argument("--sample_metadata_list", help="Path of two-column file mapping sample name to sample json "
                                                       "metadata file (for multiple samples).", default="")
    parser.add_argument("--run_metadata_list", help="Path of two-column file mapping sample name to run json file (for"
                                                    " multiple samples).", default="")
    parser.add_argument("--sample_cnv_bed_list", help="Path of two-column file mapping sample name to CNV_CALLS.bed"
                                                      " file (for multiple samples).", default="")
    args = parser.parse_args()

    output_head = ["Sample name", "AA amplicon number", "Feature ID", "Classification", "Location", "Oncogenes",
                   "All genes", "Complexity score",
                   "Captured interval length", "Feature median copy number", "Feature maximum copy number", "Filter flag", 
                   "Reference version", "Tissue of origin", "Sample type", "Feature BED file", "CNV BED file",
                   "AA PNG file", "AA PDF file", "Run metadata JSON"]

    sample_metadata_dict = defaultdict(lambda: defaultdict(lambda: "NA"))
    sample_metadata_path = defaultdict(lambda: "Not provided")
    if args.sample_metadata_list:
        with open(args.sample_metadata_list) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                sample_metadata_path[fields[0]] = fields[1]
                curr_sample_metadata = json.load(open(fields[1], 'r'))
                sample_metadata_dict[fields[0]] = curr_sample_metadata

    run_metadata_dict = defaultdict(lambda: defaultdict(lambda: "NA"))
    run_metadata_path = defaultdict(lambda: "Not provided")
    if args.run_metadata_list:
        with open(args.run_metadata_list) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                run_metadata_path[fields[0]] = fields[1]
                curr_run_metadata = json.load(open(fields[1], 'r'))
                run_metadata_dict[fields[0]] = curr_run_metadata
                # print(fields[0], curr_run_metadata)

    sample_cnv_calls_path = defaultdict(lambda: "Not provided")
    if args.sample_cnv_bed_list:
        with open(args.sample_cnv_bed_list) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                sample_cnv_calls_path[fields[0]] = fields[1]

    output_table_lines = [output_head, ]
    with open(args.input) as input_file, open(args.classification_file) as classification_file:
        classBase = args.classification_file.rsplit("_amplicon_classification_profiles.tsv")[0]
        ldir = os.path.dirname(classBase) + "/files/"
        if ldir == "/files/": ldir = "files/"

        if not os.path.exists(ldir):
            os.makedirs(ldir)

        classBedDir = classBase + "_classification_bed_files/"
        gene_file = classBase + "_gene_list.tsv"
        entropy_file = classBase + "_feature_entropy.tsv"
        basic_stats_file = classBase + "_feature_basic_properties.tsv"
        amplicon_gene_dict = read_amplicon_gene_list(gene_file)
        amplicon_complexity_dict = read_complexity_scores(entropy_file)
        basic_stats_dict = read_basic_stats(basic_stats_file)

        if args.run_metadata_file:
            metadata_dict = json.load(open(args.run_metadata_file, 'r'))
            run_metadata_dict = defaultdict(lambda: metadata_dict)
            run_metadata_path = defaultdict(lambda: os.path.abspath(args.run_metadata_file))

        if args.cnv_bed:
            if not os.path.exists(ldir + os.path.basename(args.cnv_bed)):
                shutil.copy(args.cnv_bed, ldir)

            args.cnv_bed = ldir + os.path.basename(args.cnv_bed)
            sample_cnv_calls_path = defaultdict(lambda: os.path.abspath(args.cnv_bed))

        elif args.sample_cnv_bed_list:
            for k, f in sample_cnv_calls_path.items():
                ofloc = ldir + os.path.basename(f)
                if not os.path.exists(ofloc):
                    shutil.copy(f, ldir)

                sample_cnv_calls_path[k] = os.path.abspath(ofloc)

        class_head = next(classification_file).rstrip().rsplit("\t")
        for input_line, classification_line in zip(input_file, classification_file):
            input_fields = input_line.rstrip().rsplit()
            sample_name = input_fields[0].rsplit("_amplicon")[0]
            shutil.copy(input_fields[1], ldir)
            shutil.copy(input_fields[2], ldir)
            amplicon_prefix = input_fields[1].rsplit("_cycles.txt")[0]
            if ":" not in amplicon_prefix:
                amplicon_prefix.replace("//", "/")
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

            amps_classes = []  # TODO: REFINE TO PREVENT ISSUES WITH FILTERED AMPS
            if classD["ecDNA+"] == "Positive":
                amps_classes.append(("ecDNA", int(classD["ecDNA_amplicons"])))

            if classD["BFB+"] == "Positive":
                amps_classes.append(("BFB", 1))

            elif not amps_classes and classD["amplicon_decomposition_class"] != "No amp/Invalid":
                amps_classes.append((classD["amplicon_decomposition_class"], 1))

            curr_sample_metadata = sample_metadata_dict[sample_name]
            cnv_bed_path = sample_cnv_calls_path[sample_name]
            curr_run_metadata = run_metadata_dict[sample_name]
            # print(curr_run_metadata)

            # Get the AC intervals, genes and complexity
            featureData = []
            for feature, namps in amps_classes:
                if feature == "ecDNA":
                    ecDNA_files = [x for x in os.listdir(classBedDir) if x.startswith(ampliconID) and
                                   x.rsplit("_")[-3] == "ecDNA"]

                # this is necessary since AA amps can have more than one ecDNA
                for i in range(namps):
                    if feature == "ecDNA":
                        featureBed = classBedDir + ecDNA_files[i]
                        featureID = ecDNA_files[i][:-14]

                    else:
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

                        # intervals = "|".join(interval_list)
                        intervals = str(interval_list)

                    raw_glist = amplicon_gene_dict[featureID]
                    # oncogenes = "|".join(sorted([g[0] for g in raw_glist if g[2]]))
                    oncogenes = str(sorted([g[0] for g in raw_glist if g[2]]))
                    all_genes = str(sorted([g[0] for g in raw_glist]))
                    complexity = amplicon_complexity_dict[featureID]
                    basic_stats = basic_stats_dict[featureID]

                    featureData.append([featureID, feature, intervals, oncogenes, all_genes, complexity] + basic_stats +
                                        [curr_run_metadata["ref_genome"], curr_sample_metadata["tissue_of_origin"], curr_sample_metadata["sample_type"],
                                         os.path.abspath(featureBed), cnv_bed_path])

            for ft in featureData:
                output_table_lines.append([sample_name, AA_amplicon_number] + ft + image_locs + [sample_metadata_path[sample_name],])

    tsv_ofname = classBase + "_result_table.tsv"
    # html_ofname = classBase + "_GenePatternNotebook_result_table.html"
    html_ofname = "index.html"
    json_ofname = classBase + "_result_data.json"

    with open(tsv_ofname, 'w') as outfile:
        for ll in output_table_lines:
            oline = "\t".join(ll) + "\n"
            outfile.write(oline)

    for ll in output_table_lines[1:]:
        copy_AA_files(ll, ldir)

    write_json_dict(output_table_lines, json_ofname)
    write_html_table(output_table_lines, html_ofname)
