list_of_fam="systems_done.txt"
for ref_fam in `cat ${list_of_fam}`; do
cd ${ref_fam}
echo ${ref_fam}
for i in {1..24}; do
cd anc_${i}

grep "Molecular_Weight" output.denovo_build.mol2 | awk '{print $3}' >> /gpfs/scratch/ccorbo/MW_DN_tests/MW_temp.txt
grep "Rotatable" output.denovo_build.mol2 | awk '{print $3}' >> /gpfs/scratch/ccorbo/MW_DN_tests/RB_temp.txt

cd ../
done
cd ../
done
sort -n MW_temp.txt >> Molecular_Weight_Smoothed.txt
sort -n RB_temp.txt >> Rotatable_Bond_Smoothed.txt
