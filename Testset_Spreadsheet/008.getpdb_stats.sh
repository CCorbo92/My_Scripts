cd pdb_stats
res_thresh=3.25
rfree_thresh=0.28
for sys in `cat ../clean.systems.all`
do
#insert pdb link here
wget https://www.rcsb.org/structure/${sys}
#grep the resolution
if grep -q "Resolution" ${sys}; then
  res=`cat ${sys} | grep -o -P '.Resolution.{0,25}' | grep -o -P '.>.{2,4}' | cut -c 3- | head -n1`
else
  res="n/a"
fi 

#grep the rfree
if grep -q "R-Value Free" ${sys}; then
  rfree=`cat ${sys} | grep -o -P '.R-Value Free.{0,20}' |  grep -o -P '.>.{0,10}'| cut -c 3- | head -n1`
else
  rfree="n/a"
fi

if (( $(echo "$res > $res_thresh" |bc -l) )); then
  res_flag="yes"
else
  res_flag="no"
fi

if (( $(echo "$rfree > $rfree_thresh" |bc -l) )); then
  rfree_flag="yes"
else
  rfree_flag="no"
fi

if grep -q "pharos" ${sys}; then
  pharos=`cat $sys | grep -o -P '.targets.{0,7}'| cut -c 10-| head -n1`
else
  pharos="n/a"
fi

if grep -q "uniprot" $sys; then
  uniprot=`cat $sys | grep -o -P '.uniprot/.{0,6}' | cut -c 10- | head -n1`
else
  uniprot="n/a"
fi

echo $sys " "  $res " " $rfree  " " $res_flag " " $rfree_flag " " $pharos " " $uniprot >> stats.txt 


#remove the page
rm ${sys}
done
