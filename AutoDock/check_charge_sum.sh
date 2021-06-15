list_of_fam="./../clean.systems.all"
for ref in `cat ${list_of_fam}`; do
value=`awk '{print $11}' ./${ref}/${ref}.lig.gast.pdbqt | sed '/^$/d' `
echo $value > ./${ref}/${ref}.charges.txt
echo $ref >> All_charges.txt
cd ${ref}
python ../charge_sum.py ${ref}.charges.txt >> ../All_charges.txt
cd ..
done
