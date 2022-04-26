#!/usr/bin/env python3

import os
import shutil
import sys

import pandas as pd

# takes the .input file as input (argument 1) and the amplicon_classification_profiles.tsv file as input (argument 2)

data_loc_d = {}
with open(sys.argv[1]) as infile:
    for line in infile:
        fields = line.rsplit("\t")
        sname = fields[0].rsplit("_amplicon")[0].rsplit("/")[-1]
        floc = "/".join(fields[1].rsplit("/")[:-1]) + "/"
        data_loc_d[sname] = floc


# make destinations
outloc = sys.argv[1].rsplit(".")[0] + "_images/"
if not os.path.exists(outloc):
    os.makedirs(outloc)
    for x in ["bfb/","ecdna/"]:
        for y in ["pdf","png"]:
            os.makedirs(outloc + x + y + "_softlinks/")


# read the classification results and make a softlink
classDat = pd.read_csv(sys.argv[2], sep="\t")
for ind, row in classDat.iterrows():
    sname = row["sample_name"].rsplit("_amplicon")[0].rsplit("/")[-1]
    png_fname = sname + "_" + row["amplicon_number"] + ".png"
    pdf_fname = sname + "_" + row["amplicon_number"] + ".pdf"
    data_loc = data_loc_d[sname]
    png_source_path = data_loc + png_fname
    pdf_source_path = data_loc + pdf_fname

    if row["ecDNA+"] == "Positive":
        png_dest_path = outloc + "ecdna/png_softlinks/" + png_fname
        pdf_dest_path = outloc + "ecdna/pdf_softlinks/" + png_fname
        shutil.copy(png_source_path, png_dest_path)
        # shutil.copy(pdf_source_path, pdf_dest_path)

    if row["BFB+"] == "Positive":
        png_dest_path = outloc + "bfb/png_softlinks/" + png_fname
        pdf_dest_path = outloc + "bfb/pdf_softlinks/" + png_fname
        shutil.copy(png_source_path, png_dest_path)
        # shutil.copy(pdf_source_path, pdf_dest_path)

    # elif row["BFB+"] != "Positive" and row["ecDNA+"] != "Positive":
    #
    #     if row["amplicon_classification"] == "Complex non-cyclic":
    #         shutil.copy(data_loc + "images/" + fname, outloc + "cnc/")
    #
    #     elif row["amplicon_classification"] == "Linear amplification":
    #         shutil.copy(data_loc + "images/" + fname, outloc + "linear/")
    #
    #     elif row["amplicon_classification"] == "No amp/Invalid":
    #         shutil.copy(data_loc + "images/" + fname, outloc + "invalid/")
    #
    #     else:
    #         print("BADNESS", row["amplicon_classification"], fname, outloc)