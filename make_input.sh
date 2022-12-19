#!/usr/bin/env bash

#arg 1: location to search of AA files
#arg 2: output name

touch scf.txt && rm scf.txt
touch sgf.txt && rm sgf.txt
outpre="${@: -1}"
for var in ${*%${!#}}
do
  exppath=`realpath $var`
  find $exppath -name "*_cycles.txt" | grep -v "annotated_cycles" | grep -v _classification/files/ | sort >> scf.txt
  find $exppath -name "*_graph.txt" | grep -v "features_to_graph" | grep -v "feature_to_graph" | grep -v _classification/files/ | sort >> sgf.txt
  if [ "$(wc -l < scf.txt)" -ne "$(wc -l < sgf.txt)" ]; then
    echo "ERROR: Unequal numbers of cycles and graph files found!"
    exit
  fi
done

cat scf.txt | rev | cut -f 1 -d '/' | cut -c12- | rev | sed 's/_amplicon[0-9]*$//' > san.txt
paste san.txt scf.txt sgf.txt > $outpre.input
rm san.txt scf.txt sgf.txt
