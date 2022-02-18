for i in {1..35}; do
count=`grep -w "anc_${i}" Default_Only_HMS/Orientation_Failure.txt | wc -l`
echo $i " " $count >> Default_Only_HMS/Anchor_Fail_Counts.txt
done
list_of_fam="FARMA.txt"
for ref_sys in `cat ${list_of_fam}`; do
count=`grep -w "${ref_sys}" Default_Only_HMS/Orientation_Failure.txt | wc -l`
echo $ref_sys " " $count >> Default_Only_HMS/System_Fail_Counts.txt
done
