list_of_fam="FARMA.txt"
for ref_sys in `cat ${list_of_fam}`; do
cd ${ref_sys}
echo ${ref_sys}
for i in {1..35}; do
cd anc_${i}

###Single mol2
#if ls -l | grep -q "Best_Match_Orient_Anchor.mol2" ; then
#echo "##########                                Name: ${ref_sys} anc_${i} " >> ../../Anchor_Matching_Alignment.mol2
#cat Best_Match_Orient_Anchor.mol2 >> ../../Anchor_Matching_Alignment.mol2
#echo "##########                                Name: ${ref_sys} anc_${i} " >> ../../Anchor_Matching_Alignment.mol2
#cat Selected_Anchor_reference.mol2 >> ../../Anchor_Matching_Alignment.mol2
#fi

###Mol2 for each system
if ls -l | grep -q "Best_Match_Orient_Anchor.mol2" ; then
echo "##########                                Name: ${ref_sys} anc_${i} " >> ../Anchor_Matching_Alignment.mol2
cat Best_Match_Orient_Anchor.mol2 >> ../Anchor_Matching_Alignment.mol2
echo "##########                                Name: ${ref_sys} anc_${i} " >> ../Anchor_Matching_Alignment.mol2
cat Selected_Anchor_reference.mol2 >> ../Anchor_Matching_Alignment.mol2
fi

cd ../
done
cd ../
done
