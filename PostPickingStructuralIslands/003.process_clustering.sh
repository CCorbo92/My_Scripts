# This script will best first cluster based on the flagged potential cluster pairs and the presorted list 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 11/2021
# Last Edit by: Christopher Corbo

python cluster.py PreclusterList.txt >out.txt

#How many clusters were created? This only counts clusters of greater than 1
num_clusters=`grep -v "ZINC" out.txt | tail -n1`

while read line;do
  if grep -q "${line}" out.txt ; then
  cluster="$(grep -B1 "${line}" out.txt | head -n1)"
  if echo $cluster | grep "ZINC"; then
    cluster="$(grep -B2 "${line}" out.txt | head -n1)"
  fi
  echo $line " " $cluster
  
  # if it was not clustered it still needs a cluster number and this increment starts at the number of non trivial clusters
  else
  ((num_clusters=num_clusters+1))
  echo $line " " $num_clusters 

  fi
done<"ZINC_Names_uniq.txt"
