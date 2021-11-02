# This script will insert the cluster information into header 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: --/----
# Last Edit by: Christopher Corbo
rm -r mol2_with_cluster
mkdir mol2_with_cluster
while read line;
do 
echo $line
name=`grep "Name" ./Viable/Rank_Sorted_Mol2/Viable_Sorted/xx${line} | awk '{print $3}'` 
clustnum=`grep $name clusters.txt | awk '{print $2}'`
cp ./Viable/Rank_Sorted_Mol2/Viable_Sorted/xx${line} ./mol2_with_cluster/

# 63 is the position in the header where the cluster info is inserted. If you need a different line specify it here.
sed -i "63i ##########                     Structure_Island:                   ${clustnum}" mol2_with_cluster/xx${line}
done<"viable_sorted.txt"

