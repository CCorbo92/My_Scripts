list_of_fam="smalltest10.txt"
for ref_sys in `cat ${list_of_fam}`; do
rm -r ./${ref_sys}
done
