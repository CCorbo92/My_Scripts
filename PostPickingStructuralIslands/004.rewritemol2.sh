mkdir mol2_with_cluster
while read line;
do 
name=`grep "Name" Viable/${line} | awk '{print $3}'` 
clustnum=`grep $name clusters.txt | awk '{print $2}'`
cp Viable/${line} mol2_with_cluster/

sed -i "63i ##########                     Structure_Island:                   ${clustnum}" mol2_with_cluster/${line}
done<"viable.txt"

