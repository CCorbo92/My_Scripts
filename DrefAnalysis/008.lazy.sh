list_of_fam="FARMA.txt"
for ref_sys in `cat ${list_of_fam}`; do
mkdir transfer/${ref_sys}
mv ${ref_sys}/Anchor_Matching_Alignment.mol2 transfer/${ref_sys}
done
