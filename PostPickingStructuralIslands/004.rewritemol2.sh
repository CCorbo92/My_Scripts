rm -r mol2_with_cluster
mkdir mol2_with_cluster
while read line;
do 
echo $line
name=`grep "Name" ./Viable/Rank_Sorted_Mol2/Viable_Sorted/xx${line} | awk '{print $3}'` 
clustnum=`grep $name clusters.txt | awk '{print $2}'`
cp ./Viable/Rank_Sorted_Mol2/Viable_Sorted/xx${line} ./mol2_with_cluster/

sed -i "63i ##########                     Structure_Island:                   ${clustnum}" mol2_with_cluster/xx${line}
done<"viable_sorted.txt"

