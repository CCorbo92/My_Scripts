# This script fills in the crossdock cvs for systems which did not have a family and fills them in as "-" to help with final alignments
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 08/2021
# Last Edit by: Christopher Corbo

#Lines in CrossDock.csv minus one
testset_size="742"

#sorting is necessary for diff
sort -n CrossDock.csv| tail -n  ${testset_size} | awk -F',' '{print $1}' | rev | cut -c 1- | rev > CrossDock_Sort.csv 

#Ignore trailing white space and diff to find which systems from the testset are not in crossdocking
diff -Z CrossDock_Sort.csv clean.systems.all | grep ">" | cut -c 3- > Not_In_Crossdock.txt

#Append a space holder for the systems in Not_In_Crossdock.txt
while read line; do
echo ${line} ",-,-,-" >> Not_In_Crossdock.csv
done <"Not_In_Crossdock.txt"

cat Not_In_Crossdock.csv >> CrossDock.csv

sort -n CrossDock.csv > CrossDock_All_sys.csv
