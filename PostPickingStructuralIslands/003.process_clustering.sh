
num_clusters=$1
python hms_clusterv2.py PreclusterList.txt >out.txt
while read line;do
  #echo $line 
  if grep -q "${line}" out.txt ; then
  cluster="$(grep -B1 "${line}" out.txt | head -n1)"
  if echo $cluster | grep "ZINC"; then
    cluster="$(grep -B2 "${line}" out.txt | head -n1)"
  fi
  echo $line " " $cluster

  else
  ((num_clusters=num_clusters+1))
  echo $line " " $num_clusters 

  fi
done<"ZINC_Names_uniq.txt"
