# This script will check for any cases in growth trees with dref and std score strings if the full ref did better than the dref assignment in any cases 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: --/----
# Last Edit by: Christopher Corbo
cd 2GQG

for i in {1..22}; do
cd anc_$i

echo $i

ls | grep "output.growth_tree_" >> list_gtree.txt
for file in `cat list_gtree.txt`; do
dref=`grep "Score_String_DRef" $file | awk -F',' '{print $NF}' | head -n1`
sref=`grep "Score_String_Std" $file | awk -F',' '{print $NF}' | head -n1`

echo $dref " " $sref >> ../ref_comp.txt

done

python ../012.compscript.py

cd ..
done
