#!/usr/bin/env bash

#arg 1: location to search of AA files
#arg 2: output name

# first, find all the AA graph and cycles files from AA runs.
rm -f scf.txt
rm -f sgf.txt
outpre="${@: -1}"
for var in ${*%${!#}}
do
  exppath=`realpath $var`
  find $exppath -name "*_cycles.txt" | grep -v "annotated_cycles" | grep -v "/\._" | grep -v "BPG_converted" | grep -v _classification/files/ | sort >> scf.txt
  find $exppath -name "*_graph.txt" | grep -v "features_to_graph" | grep -v "/\._" | grep -v "feature_to_graph" | grep -v _classification/files/ | sort >> sgf.txt
  if [ "$(wc -l < scf.txt)" -ne "$(wc -l < sgf.txt)" ]; then
    echo "ERROR: Unequal numbers of cycles and graph files found! AA may not have completed correctly."
    exit 1
  fi
done

cat scf.txt | rev | cut -f 1 -d '/' | cut -c12- | rev | sed 's/_amplicon[0-9]*$//' > san.txt
paste san.txt scf.txt sgf.txt > $outpre.input

# now find the _summary.txt files so that samples are included which did not have corresponding AA outputs.
rm -f ssf.txt
for var in ${*%${!#}}
do
  exppath=`realpath $var`
  find $exppath -name "*_summary.txt" | grep -v "/\._" | grep -v _classification/files/ | sort >> ssf.txt
done

cat ssf.txt | rev | cut -f 1 -d '/' | cut -c13- | rev | sed 's/_summary.txt$//' > ssn.txt
paste ssn.txt ssf.txt > ${outpre}_summary_map.txt

rm -f san.txt scf.txt sgf.txt ssf.txt ssn.txt

