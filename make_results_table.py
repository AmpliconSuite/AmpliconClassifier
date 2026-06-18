#!/usr/bin/env python3

import argparse
from collections import defaultdict
import json
import logging
import os
import shutil
import sys


NON_FEATURE_CLASSES = {"No amp/Invalid", "No-FSCNA", "Invalid"}
FALLBACK_FEATURE_CLASSES = {"Linear", "Complex-non-cyclic"}
MECHANISM_DECOMPOSITION_CLASSES = {"Virus"}
OUTPUT_HEAD = ["Sample name", "AA amplicon number", "Feature ID", "Classification",
               "FAN probability", "Location", "Oncogenes", "All genes", "NCBI Gene IDs",
               "Complexity score", "ecDNA context", "Captured interval length", "Feature median copy number",
               "Feature maximum copy number", "Filter flag", "Reference version", "Tissue of origin",
               "Sample type", "Feature BED file", "CNV BED file", "AS-p version", "AA version", "AC version",
               "AA PNG file", "AA PDF file", "AA summary file", "Run metadata JSON", "Sample metadata JSON"]
FILE_COLUMNS = ["AA PNG file", "AA PDF file", "AA summary file", "Run metadata JSON", "Sample metadata JSON"]


def parse_bool(value):
    if isinstance(value, bool):
        return value

    sval = str(value).strip().lower()
    if sval in {"true", "t", "1", "yes", "y"}:
        return True
    if sval in {"false", "f", "0", "no", "n", "", "na", "none"}:
        return False

    raise ValueError("Cannot parse boolean value: {!r}".format(value))


def feature_sort_key(fname):
    base = fname.rsplit("_intervals.bed", 1)[0]
    suffix = base.rsplit("_", 1)[-1]
    try:
        return int(suffix)
    except ValueError:
        return suffix


def copy_if_needed(src, dest_dir):
    if src in {"Not found", "Not provided", "NA"} or src.endswith("Not found") or src.endswith("Not provided"):
        return src

    dest = os.path.join(dest_dir, os.path.basename(src))
    if not os.path.exists(dest):
        shutil.copy(src, dest_dir)

    return dest


def copy_result_files(row, output_dir):
    for col in FILE_COLUMNS:
        row[col] = copy_if_needed(row[col], output_dir)


def row_values(row):
    return [str(row.get(col, "NA")) for col in OUTPUT_HEAD]


def write_result_outputs(rows, tsv_ofname, json_ofname, files_dir):
    with open(tsv_ofname, 'w') as outfile:
        outfile.write("\t".join(OUTPUT_HEAD) + "\n")
        for row in rows:
            outfile.write("\t".join(row_values(row)) + "\n")

    for row in rows:
        copy_result_files(row, files_dir)

    with open(json_ofname, 'w') as outfile:
        json.dump([{col: row.get(col, "NA") for col in OUTPUT_HEAD} for row in rows],
                  outfile, sort_keys=True, indent=2)


def read_amplicon_gene_list(gene_file):
    amplicon_gene_dict = defaultdict(list)
    with open(gene_file) as infile:
        h = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            fd = dict(zip(h, fields))
            featureID = "_".join(fields[:3])
            if 'ncbi_id' in fd:
                amplicon_gene_dict[featureID].append(
                    (fd['gene'], fd['gene_cn'], parse_bool(fd['is_canonical_oncogene']), fd['ncbi_id']))
            else:
                amplicon_gene_dict[featureID].append(
                    (fd['gene'], fd['gene_cn'], parse_bool(fd['is_canonical_oncogene']), "NA"))

    return amplicon_gene_dict


def read_complexity_scores(complexity_file):
    amplicon_complexity_dict = defaultdict(lambda: "NA")
    with open(complexity_file) as infile:
        h = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            fd = dict(zip(h, fields))
            featureID = "_".join(fields[:3])
            amplicon_complexity_dict[featureID] = fd["feature_complexity"]

    return amplicon_complexity_dict


def read_context(context_file):
    amplicon_context_dict = defaultdict(lambda: "NA")
    if not os.path.exists(context_file):
        sys.stderr.write("Could not locate ecDNA context calls file. Classification results may be outdated.")

    else:
        with open(context_file) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                featureID = fields[0]
                amplicon_context_dict[featureID] = fields[1]

    return amplicon_context_dict


def read_basic_stats(basic_stats_file):
    basic_stats_dict = defaultdict(lambda: ["NA", "NA", "NA", "NA"])
    with open(basic_stats_file) as infile:
        h = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            featureID = fields[0]
            basic_stats_dict[featureID] = fields[1:]

    return basic_stats_dict


def read_fan_calls(fan_file):
    fan_dict = defaultdict(lambda: {"call": "NA", "probability": "NA"})
    if not os.path.exists(fan_file):
        return fan_dict

    with open(fan_file) as infile:
        try:
            h = next(infile).rstrip().rsplit("\t")
        except StopIteration:
            return fan_dict
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            fd = dict(zip(h, fields))
            amplicon_id = "_".join([fd["sample_name"], fd["amplicon_number"]])
            fan_dict[amplicon_id] = {
                "call": fd.get("fan_call", "NA"),
                "probability": fd.get("fan_probability", "NA"),
            }

    return fan_dict


def read_summary_list(summ_map_file):
    sumf_dict = {}
    if not summ_map_file:
        return sumf_dict

    if not os.path.exists(summ_map_file):
        sys.stderr.write("Warning: summary map {} not found. Continuing without AA summary files.\n".format(
            summ_map_file))
        return sumf_dict

    # key is the sample name plus the base path of the summary file.
    with open(summ_map_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            floc = os.path.dirname(fields[1])
            sumf_dict[(fields[0], floc)] = fields[1]

    return sumf_dict


def cycles_graph_amp_lookup(input_file_path):
    lookup = {}
    with open(input_file_path) as input_file:
        for line_num, line in enumerate(input_file, 1):
            if not line.strip():
                continue
            fields = line.rstrip().rsplit()
            if len(fields) < 3:
                raise ValueError("Malformed input line {} in {}: {}".format(line_num, input_file_path, line.rstrip()))
            a_id = fields[1].rsplit("/")[-1].rsplit("_cycles.txt")[0]
            if a_id in lookup:
                raise ValueError("Duplicate amplicon ID {} in {}".format(a_id, input_file_path))
            lookup[a_id] = (fields[1], fields[2])

    return lookup


def get_version_info(run_meta_dict):
    asp, aa, ac = "NA", "NA", "NA"
    if "AC_version" in run_meta_dict and run_meta_dict["AC_version"]:
        ac = run_meta_dict["AC_version"]

    if "PAA_version" in run_meta_dict and run_meta_dict["PAA_version"]:
        asp = run_meta_dict["PAA_version"]

    if "AmpliconSuite-pipeline_version" in run_meta_dict and run_meta_dict["AmpliconSuite-pipeline_version"]:
        asp = run_meta_dict["AmpliconSuite-pipeline_version"]

    if "AA_version" in run_meta_dict and run_meta_dict["AA_version"]:
        aa = run_meta_dict["AA_version"].rsplit()[-1]

    return asp, aa, ac


def read_bed_intervals(feature_bed):
    if not os.path.exists(feature_bed):
        sys.stderr.write("Warning: interval file " + feature_bed + " not found!\n")
        return "Interval file not found"

    interval_list = []
    with open(feature_bed) as bedfile:
        for line in bedfile:
            fields = line.rstrip().rsplit("\t")
            if len(fields) >= 3:
                interval_list.append(fields[0] + ":" + fields[1] + "-" + fields[2])

    return str(interval_list)


def discover_feature_specs(classD, amplicon_id, class_bed_dir):
    feature_specs = []
    if classD["ecDNA+"] == "Positive":
        ecDNA_files = sorted(
            [x for x in os.listdir(class_bed_dir)
             if x.startswith(amplicon_id + "_") and x.rsplit("_")[-3] == "ecDNA"],
            key=feature_sort_key
        )
        for ind in range(int(classD["ecDNA_amplicons"])):
            if ind >= len(ecDNA_files):
                sys.stderr.write("Warning: expected ecDNA interval file for {} ecDNA_{} not found!\n".format(
                    amplicon_id, ind + 1))
                feature_id = "_".join([amplicon_id, "ecDNA", str(ind + 1)])
                feature_bed = os.path.join(class_bed_dir, feature_id + "_intervals.bed")
            else:
                feature_bed = os.path.join(class_bed_dir, ecDNA_files[ind])
                feature_id = ecDNA_files[ind][:-14]
            feature_specs.append((feature_id, "ecDNA", feature_bed))

    if classD["BFB+"] == "Positive":
        feature_id = "_".join([amplicon_id, "BFB", "1"])
        feature_specs.append((feature_id, "BFB", os.path.join(class_bed_dir, feature_id + "_intervals.bed")))

    if classD.get("FAN+", "None detected") == "Positive":
        feature_id = "_".join([amplicon_id, "FAN", "1"])
        feature_specs.append((feature_id, "FAN", os.path.join(class_bed_dir, feature_id + "_intervals.bed")))

    feature = classD["amplicon_decomposition_class"]
    if feature in MECHANISM_DECOMPOSITION_CLASSES:
        feature_id = "_".join([amplicon_id, feature, "1"])
        feature_specs.append((feature_id, feature, os.path.join(class_bed_dir, feature_id + "_intervals.bed")))
    elif not feature_specs and feature in FALLBACK_FEATURE_CLASSES:
        feature_id = "_".join([amplicon_id, feature, "1"])
        feature_specs.append((feature_id, feature, os.path.join(class_bed_dir, feature_id + "_intervals.bed")))
    elif not feature_specs and feature == "Cyclic":
        logging.warning(
            "Amplicon %s has decomposition class Cyclic but no ecDNA, BFB, FAN, or Virus mechanism; "
            "not reporting a Cyclic feature row.",
            amplicon_id,
        )

    return feature_specs


def build_feature_row(sample_name, amplicon_number, feature_id, feature, feature_bed, chromo_prob, intervals,
                      gene_dict, complexity_dict, context_dict, basic_stats_dict, run_metadata, sample_metadata,
                      cnv_bed_path, versions, image_locs, summary_files):
    sorted_glist = sorted(gene_dict[feature_id])
    asp_version, aa_version, ac_version = versions
    sumf, run_metadata_file, sample_metadata_file = summary_files
    basic_stats = basic_stats_dict[feature_id]
    return {
        "Sample name": sample_name,
        "AA amplicon number": amplicon_number,
        "Feature ID": feature_id,
        "Classification": feature,
        "FAN probability": chromo_prob,
        "Location": intervals,
        "Oncogenes": str([g[0] for g in sorted_glist if g[2]]),
        "All genes": str([g[0] for g in sorted_glist]),
        "NCBI Gene IDs": str([g[3] for g in sorted_glist]),
        "Complexity score": complexity_dict[feature_id],
        "ecDNA context": context_dict[feature_id],
        "Captured interval length": basic_stats[0],
        "Feature median copy number": basic_stats[1],
        "Feature maximum copy number": basic_stats[2],
        "Filter flag": basic_stats[3],
        "Reference version": run_metadata["ref_genome"],
        "Tissue of origin": sample_metadata["tissue_of_origin"],
        "Sample type": sample_metadata["sample_type"],
        "Feature BED file": os.path.abspath(feature_bed),
        "CNV BED file": cnv_bed_path,
        "AS-p version": asp_version,
        "AA version": aa_version,
        "AC version": ac_version,
        "AA PNG file": image_locs[0],
        "AA PDF file": image_locs[1],
        "AA summary file": sumf,
        "Run metadata JSON": run_metadata_file,
        "Sample metadata JSON": sample_metadata_file,
    }


def build_no_feature_row(sample_name, sumf, sample_metadata, run_metadata, cnv_bed_path, versions, summary_files):
    asp_version, aa_version, ac_version = versions
    _, run_metadata_file, sample_metadata_file = summary_files
    feature_id = sample_name + "_NA"
    basic_stats = ["NA", "NA", "NA", "NA"]
    return {
        "Sample name": sample_name,
        "AA amplicon number": "NA",
        "Feature ID": feature_id,
        "Classification": "NA",
        "FAN probability": "NA",
        "Location": "[]",
        "Oncogenes": "[]",
        "All genes": "[]",
        "NCBI Gene IDs": "[]",
        "Complexity score": "NA",
        "ecDNA context": "NA",
        "Captured interval length": basic_stats[0],
        "Feature median copy number": basic_stats[1],
        "Feature maximum copy number": basic_stats[2],
        "Filter flag": basic_stats[3],
        "Reference version": run_metadata["ref_genome"],
        "Tissue of origin": sample_metadata["tissue_of_origin"],
        "Sample type": sample_metadata["sample_type"],
        "Feature BED file": os.path.abspath("NA"),
        "CNV BED file": cnv_bed_path,
        "AS-p version": asp_version,
        "AA version": aa_version,
        "AC version": ac_version,
        "AA PNG file": "NA",
        "AA PDF file": "NA",
        "AA summary file": sumf,
        "Run metadata JSON": run_metadata_file,
        "Sample metadata JSON": sample_metadata_file,
    }


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


def make_results_table(input_file, classification_file, summary_map=None, sample_metadata_file="",
                       run_metadata_file="", cnv_bed="", sample_metadata_list="", run_metadata_list="",
                       sample_cnv_bed_list="", ref=None):
    """
    Main function to create results table from AC outputs.

    Args:
        input_file: Path to .input file produced by make_input.sh
        classification_file: Path of amplicon_classification_profiles.tsv file
        summary_map: Path to the _summary_map.txt file (optional)
        sample_metadata_file: Path of sample metadata json (optional)
        run_metadata_file: Path of run metadata json (optional)
        cnv_bed: Path of the CNV_CALLS.bed file (optional)
        sample_metadata_list: Path of file mapping sample name to metadata json (optional)
        run_metadata_list: Path of file mapping sample name to run json (optional)
        sample_cnv_bed_list: Path of file mapping sample name to CNV bed (optional)
        ref: Reference genome name (optional if metadata files provided)
    """
    if not run_metadata_file and not run_metadata_list and not ref:
        sys.stderr.write("One of the following must be provided: --ref | --run_metadata_list | --run_metadata_file\n")
        sys.exit(1)

    classBase = classification_file.rsplit("_amplicon_classification_profiles.tsv")[0]
    ldir = os.path.dirname(classBase) + "/files/"
    if ldir == "/files/": ldir = "files/"
    if not os.path.exists(ldir): os.makedirs(ldir)

    if not summary_map:
        summary_map = input_file.rsplit(".input", 1)[0] + "_summary_map.txt"

    sumf_used = set()
    sumf_dict = read_summary_list(summary_map)
    sample_metadata_dict = defaultdict(lambda: defaultdict(lambda: "NA"))
    sample_metadata_path = defaultdict(lambda: "Not provided")
    if sample_metadata_list:
        with open(sample_metadata_list) as infile:
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
    if run_metadata_list:
        with open(run_metadata_list) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                crmf = os.path.abspath(fields[1])
                run_metadata_path[fields[0]] = os.path.abspath(crmf)
                curr_run_metadata = json.load(open(crmf, 'r'))
                run_metadata_dict[fields[0]] = curr_run_metadata
                if not os.path.exists(ldir + os.path.basename(crmf)):
                    shutil.copy(crmf, ldir)

    sample_cnv_calls_path = defaultdict(lambda: "Not provided")
    if sample_cnv_bed_list:
        with open(sample_cnv_bed_list) as infile:
            for line in infile:
                fields = line.rstrip().rsplit()
                sample_cnv_calls_path[fields[0]] = fields[1]

    result_rows = []
    cyc_graph_lookup_dct = cycles_graph_amp_lookup(input_file)
    with open(classification_file) as classification_file_handle:
        classBedDir = classBase + "_classification_bed_files/"
        gene_file = classBase + "_gene_list.tsv"
        complexity_file = classBase + "_feature_complexity.tsv"
        basic_stats_file = classBase + "_feature_basic_properties.tsv"
        context_file = classBase + "_ecDNA_context_calls.tsv"
        fan_file = classBase + "_fan_calls.tsv"
        amplicon_gene_dict = read_amplicon_gene_list(gene_file)
        amplicon_complexity_dict = read_complexity_scores(complexity_file)
        basic_stats_dict = read_basic_stats(basic_stats_file)
        context_dict = read_context(context_file)
        fan_dict = read_fan_calls(fan_file)

        if sample_metadata_file:
            init_sample_metadata = json.load(open(sample_metadata_file, 'r'))
            sample_metadata_dict = defaultdict(lambda: init_sample_metadata)
            sample_metadata_path = defaultdict(lambda: os.path.abspath(sample_metadata_file))
            if not os.path.exists(ldir + os.path.basename(sample_metadata_file)):
                shutil.copy(sample_metadata_file, ldir)

        if run_metadata_file:
            init_run_metadata = json.load(open(run_metadata_file, 'r'))
            run_metadata_dict = defaultdict(lambda: init_run_metadata)
            run_metadata_path = defaultdict(lambda: os.path.abspath(run_metadata_file))
            if not os.path.exists(ldir + os.path.basename(run_metadata_file)):
                shutil.copy(run_metadata_file, ldir)

        if cnv_bed:
            if not os.path.exists(ldir + os.path.basename(cnv_bed)):
                shutil.copy(cnv_bed, ldir)

            cnv_bed = ldir + os.path.basename(cnv_bed)
            sample_cnv_calls_path = defaultdict(lambda: os.path.abspath(cnv_bed))

        elif sample_cnv_bed_list:
            for k, f in sample_cnv_calls_path.items():
                ofloc = ldir + os.path.basename(f)
                if not os.path.exists(ofloc):
                    shutil.copy(f, ldir)

                sample_cnv_calls_path[k] = os.path.abspath(ofloc)

        class_head = next(classification_file_handle).rstrip().rsplit("\t")
        for classification_line in classification_file_handle:
            classD = dict(zip(class_head, classification_line.rstrip().rsplit("\t")))
            sample_name = classD['sample_name']
            ampliconID = "_".join([classD["sample_name"], classD["amplicon_number"]])
            if ampliconID not in cyc_graph_lookup_dct:
                raise ValueError("Classification amplicon {} not found in input {}".format(ampliconID, input_file))
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
                sys.stderr.write("Amplicon names in " + input_file + " do not match "
                                 + classification_file + "\n")
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

            curr_sample_metadata = sample_metadata_dict[sample_name]
            cnv_bed_path = sample_cnv_calls_path[sample_name]
            curr_run_metadata = run_metadata_dict[sample_name]
            if curr_run_metadata['ref_genome'] == "NA" and ref:
                curr_run_metadata['ref_genome'] = ref
            versions = get_version_info(curr_run_metadata)
            summary_files = (sumf, run_metadata_path[sample_name], sample_metadata_path[sample_name])

            feature_specs = discover_feature_specs(classD, ampliconID, classBedDir)
            if feature_specs:
                sumf_used.add((sample_name, cfile_dir))
            chromo_prob = fan_dict[ampliconID]["probability"]
            for featureID, feature, featureBed in feature_specs:
                result_rows.append(build_feature_row(
                    sample_name, AA_amplicon_number, featureID, feature, featureBed, chromo_prob,
                    read_bed_intervals(featureBed), amplicon_gene_dict, amplicon_complexity_dict, context_dict,
                    basic_stats_dict, curr_run_metadata, curr_sample_metadata, cnv_bed_path, versions, image_locs,
                    summary_files
                ))

        for k in set(sumf_dict.keys()) - sumf_used:
            print(k[0] + " had no identifiable focal amplifications in the AA amplicons")
            sumf = sumf_dict[k]
            sample_name = k[0]
            curr_sample_metadata = sample_metadata_dict[sample_name]
            cnv_bed_path = sample_cnv_calls_path[sample_name]
            curr_run_metadata = run_metadata_dict[sample_name]
            if curr_run_metadata['ref_genome'] == "NA" and ref:
                curr_run_metadata["ref_genome"] = ref

            versions = get_version_info(curr_run_metadata)
            summary_files = (sumf, run_metadata_path[sample_name], sample_metadata_path[sample_name])
            result_rows.append(build_no_feature_row(
                sample_name, sumf, curr_sample_metadata, curr_run_metadata, cnv_bed_path, versions, summary_files
            ))

    tsv_ofname = classBase + "_result_table.tsv"
    # html_ofname = "index.html"
    json_ofname = classBase + "_result_data.json"

    write_result_outputs(result_rows, tsv_ofname, json_ofname, ldir)
    print("Finished creating summary tables for " + str(len(result_rows)) + " total entries")


if __name__ == "__main__":
    # The input file must be the source of the classification file calls.
    parser = argparse.ArgumentParser(description="Organize AC results into a table")
    parser.add_argument("-i", "--input", help="Path to .input file produced by make_input.sh. Each line formatted as: "
                                              "sample_name cycles.txt graph.txt (required)", required=True)
    parser.add_argument("--classification_file", help="Path of amplicon_classification_profiles.tsv file (required)",
                        required=True)
    parser.add_argument("--summary_map", help="Path to the _summary_map.txt file produced by make_input.sh")
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

    # Call the main function
    make_results_table(
        input_file=args.input,
        classification_file=args.classification_file,
        summary_map=args.summary_map,
        sample_metadata_file=args.sample_metadata_file,
        run_metadata_file=args.run_metadata_file,
        cnv_bed=args.cnv_bed,
        sample_metadata_list=args.sample_metadata_list,
        run_metadata_list=args.run_metadata_list,
        sample_cnv_bed_list=args.sample_cnv_bed_list,
        ref=args.ref
    )
