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
    basic_stats_dict = defaultdict(lambda: ["NA", "NA", "NA", "NA"])
    with open(basic_stats_file) as infile:
        h = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            featureID = fields[0]
            basic_stats_dict[featureID] = fields[1:]

    return basic_stats_dict


def read_summary_list(summ_map_file):
    sumf_dict = {}
    if not summ_map_file:
        return sumf_dict

    # key is the sample name plus the base path of the summary file.
    with open(summ_map_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            floc = os.path.dirname(fields[1])
            sumf_dict[(fields[0], floc)] = fields[1]

    return sumf_dict


def cycles_graph_amp_lookup(input_lines):
    lookup = {}
    with open(args.input) as input_file:
        for line in input_file:
            fields = line.rstrip().rsplit()
            a_id = fields[1].rsplit("/")[-1].rsplit("_cycles.txt")[0]
            lookup[a_id] = (fields[1], fields[2])

    return lookup


def copy_AA_files(ll, ldir):
    for i in range(-5, 0):
        s = ll[i]
        if not s.endswith("Not found") and not s.endswith("Not provided") and not s == "NA":
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
        json.dump(dlist, outfile, sort_keys=True, indent=2)

    pass


if __name__ == "__main__":
    # The input file must be the source of the classification file calls.
    parser = argparse.ArgumentParser(description="Organize AC results into a table")
    parser.add_argument("-i", "--input", help="Path to .input file produced by make_input.sh. Each line formatted as: "
                        "sample_name cycles.txt graph.txt.", required=True)
    parser.add_argument("--classification_file", help="Path of amplicon_classification_profiles.tsv file",
                        required=True)
    parser.add_argument("--summary_map", help="Path to the _summary_map.txt file produced by make_input.sh",
                        required=True)
    parser.add_argument("--sample_metadata_file", help="Path of sample metadata, [sample]_sample_metadata.json file"
                                                       " (for building table with a single sample).", default="")
    parser.add_argument("--run_metadata_file", help="Path of run metadata, [sample]_run_metadata.json file (for "
                                                    "building table with a single sample).", default="")
    parser.add_argument("--cnv_bed", help="Path of the CNV_CALLS.bed file used for this run (for single sample).",
                        default="")
    parser.add_argument("--sample_metadata_list", help="Path of two-column file mapping sample name to sample json "
                                                       "metadata file (for multiple samples).", default="")
    parser.add_argument("--run_metadata_list", help="Path of two-column file mapping sample name to run json file (for"
                                                    " multiple samples).", default="")
    parser.add_argument("--sample_cnv_bed_list", help="Path of two-column file mapping sample name to CNV_CALLS.bed"
                                                      " file (for multiple samples).", default="")
    parser.add_argument("--ref", help="Reference genome name used for alignment, one of hg19, GRCh37, or GRCh38. "
                                      "Required if --run_metadata_file or --run_metadata_list are not set",
                        choices=["hg19", "GRCh37", "GRCh38", "GRCh38_viral", "mm10"])
    args = parser.parse_args()

    if not args.run_metadata_file and not args.run_metadata_list and not args.ref:
        sys.stderr.write("One of the following must be provided: --ref | --run_metadata_list | --run_medata_file\n")

    output_head = ["Sample name", "AA amplicon number", "Feature ID", "Classification", "Location", "Oncogenes",
                   "All genes", "Complexity score", "Captured interval length", "Feature median copy number",
                   "Feature maximum copy number", "Filter flag", "Reference version", "Tissue of origin",
                   "Sample type", "Feature BED file", "CNV BED file", "AA PNG file", "AA PDF file", "AA summary file",
                   "Run metadata JSON", "Sample metadata JSON"]

    sumf_used = set()
    sumf_dict = read_summary_list(args.summary_map)

    classBase = args.classification_file.rsplit("_amplicon_classification_profiles.tsv")[0]
    ldir = os.path.dirname(classBase) + "/files/"
    if ldir == "/files/": ldir = "files/"
    if not os.path.exists(ldir): os.makedirs(ldir)

    sample_metadata_dict = defaultdict(lambda: defaultdict(lambda: "NA"))
    sample_metadata_path = defaultdict(lambda: "Not provided")
    if args.sample_metadata_list:
        with open(args.sample_metadata_list) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                csmf = os.path.abspath(fields[1])
                sample_metadata_path[fields[0]] = os.path.abspath(csmf)
                curr_sample_metadata = json.load(open(csmf, 'r'))
                sample_metadata_dict[fields[0]] = curr_sample_metadata
                if not os.path.exists(ldir + os.path.basename(csmf)):
                    shutil.copy(csmf, ldir)

    run_metadata_dict = defaultdict(lambda: defaultdict(lambda: "NA"))
    run_metadata_path = defaultdict(lambda: "Not provided")
    if args.run_metadata_list:
        with open(args.run_metadata_list) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                crmf = os.path.abspath(fields[1])
                run_metadata_path[fields[0]] = os.path.abspath(crmf)
                curr_run_metadata = json.load(open(crmf, 'r'))
                run_metadata_dict[fields[0]] = curr_run_metadata
                if not os.path.exists(ldir + os.path.basename(crmf)):
                    shutil.copy(crmf, ldir)

    sample_cnv_calls_path = defaultdict(lambda: "Not provided")
    if args.sample_cnv_bed_list:
        with open(args.sample_cnv_bed_list) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                sample_cnv_calls_path[fields[0]] = fields[1]

    output_table_lines = [output_head, ]
    cyc_graph_lookup_dct = cycles_graph_amp_lookup(args.input)
    with open(args.classification_file) as classification_file:
        classBedDir = classBase + "_classification_bed_files/"
        gene_file = classBase + "_gene_list.tsv"
        entropy_file = classBase + "_feature_entropy.tsv"
        basic_stats_file = classBase + "_feature_basic_properties.tsv"
        amplicon_gene_dict = read_amplicon_gene_list(gene_file)
        amplicon_complexity_dict = read_complexity_scores(entropy_file)
        basic_stats_dict = read_basic_stats(basic_stats_file)

        if args.sample_metadata_file:
            init_sample_metadata = json.load(open(args.sample_metadata_file, 'r'))
            sample_metadata_dict = defaultdict(lambda: init_sample_metadata)
            sample_metadata_path = defaultdict(lambda: os.path.abspath(args.sample_metadata_file))
            if not os.path.exists(ldir + os.path.basename(args.sample_metadata_file)):
                shutil.copy(args.sample_metadata_file, ldir)

        if args.run_metadata_file:
            init_run_metadata = json.load(open(args.run_metadata_file, 'r'))
            run_metadata_dict = defaultdict(lambda: init_run_metadata)
            run_metadata_path = defaultdict(lambda: os.path.abspath(args.run_metadata_file))
            if not os.path.exists(ldir + os.path.basename(args.run_metadata_file)):
                shutil.copy(args.run_metadata_file, ldir)

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
        for classification_line in classification_file:
            classD = dict(zip(class_head, classification_line.rstrip().rsplit("\t")))
            sample_name = classD['sample_name']
            ampliconID = "_".join([classD["sample_name"], classD["amplicon_number"]])
            cycles_file, graph_file = cyc_graph_lookup_dct[ampliconID]
            shutil.copy(cycles_file, ldir)
            shutil.copy(graph_file, ldir)

            # what is the directory of the cycles file?
            cfile_dir = os.path.dirname(cycles_file)
            if (sample_name, cfile_dir) in sumf_dict:
                sumf = sumf_dict[(sample_name, cfile_dir)]
            else:
                sumf = "Not provided"

            amplicon_prefix = cycles_file.rsplit("_cycles.txt")[0]
            if ":" not in amplicon_prefix:
                amplicon_prefix.replace("//", "/")

            AA_amplicon_number = classD["amplicon_number"].lstrip("amplicon")

            # TODO: REVISE THIS ERROR AFTER NEW MATCHING OF INPUTS TO CLASSIFICATIONS
            if sample_name + "_amplicon" + AA_amplicon_number != ampliconID:
                sys.stderr.write(sample_name + "_amplicon" + AA_amplicon_number + " | " + ampliconID + "\n")
                sys.stderr.write("Amplicon names in " + args.input + " do not with "
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
                sumf_used.add((sample_name, cfile_dir))
                amps_classes.append(("ecDNA", int(classD["ecDNA_amplicons"])))

            if classD["BFB+"] == "Positive":
                sumf_used.add((sample_name, cfile_dir))
                amps_classes.append(("BFB", 1))

            elif not amps_classes and classD["amplicon_decomposition_class"] != "No amp/Invalid":
                sumf_used.add((sample_name, cfile_dir))
                amps_classes.append((classD["amplicon_decomposition_class"], 1))

            curr_sample_metadata = sample_metadata_dict[sample_name]
            cnv_bed_path = sample_cnv_calls_path[sample_name]
            curr_run_metadata = run_metadata_dict[sample_name]
            if curr_run_metadata['ref_genome'] == "NA" and args.ref:
                curr_run_metadata['ref_genome'] = args.ref
            # print(curr_run_metadata)

            # Get the AC intervals, genes and complexity
            featureData = []
            for feature, namps in amps_classes:
                if feature == "ecDNA":
                    ecDNA_files = [x for x in os.listdir(classBedDir) if x.startswith(ampliconID + "_") and
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
                output_table_lines.append(
                    [sample_name, AA_amplicon_number] + ft + image_locs +
                    [sumf, run_metadata_path[sample_name], sample_metadata_path[sample_name]])

        for k in set(sumf_dict.keys()) - sumf_used:
            print(k[0] + " had no identifiable focal amplifications in the AA amplicons")
            sumf = sumf_dict[k]
            sample_name = k[0]
            AA_amplicon_number = "NA"
            feature = "NA"
            featureID = sample_name + "_" + feature
            intervals = "[]"
            oncogenes = "[]"
            all_genes = "[]"
            complexity = "NA"
            basic_stats = basic_stats_dict[featureID]
            featureBed = "NA"
            curr_sample_metadata = sample_metadata_dict[sample_name]
            cnv_bed_path = sample_cnv_calls_path[sample_name]
            curr_run_metadata = run_metadata_dict[sample_name]
            if curr_run_metadata['ref_genome'] == "NA" and args.ref:
                curr_run_metadata["ref_genome"] = args.ref

            fdl = [featureID, feature, intervals, oncogenes, all_genes, complexity] + basic_stats + \
                  [curr_run_metadata["ref_genome"], curr_sample_metadata["tissue_of_origin"],
                   curr_sample_metadata["sample_type"], os.path.abspath(featureBed), cnv_bed_path]

            image_locs = ["NA", "NA"]
            output_table_lines.append(
                [sample_name, AA_amplicon_number] + fdl + image_locs +
                [sumf, run_metadata_path[sample_name], sample_metadata_path[sample_name], ])


    tsv_ofname = classBase + "_result_table.tsv"
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
    print("Finished creating summary tables for " + str(len(output_table_lines[1:])) + " total entries")
