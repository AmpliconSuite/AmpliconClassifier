#!/usr/bin/env bash

#arg 1: location to search of AA files
#arg 2: output name

exppath=`realpath $1`
find $exppath -name "*_cycles.txt" | grep -v "annotated_cycles" | sort > scf.txt
find $exppath -name "*_graph.txt" | grep -v "features_to_graph" | grep -v "feature_to_graph" | sort > sgf.txt
if [ "$(wc -l < scf.txt)" -ne "$(wc -l < sgf.txt)" ]; then
  echo "ERROR: Unequal numbers of cycles and graph files found!"
  exit
fi
cat scf.txt | rev | cut -d '/' -f 1 | cut -c12- | rev | sed 's/_amplicon[0-9]*$//' > san.txt
paste san.txt scf.txt sgf.txt > $2.input
rm san.txt scf.txt sgf.txt
