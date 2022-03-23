cd 2GQG

grep -A1 "USER_CHARGES" 2GQG_rescored_fullref_ranked.mol2 | grep "\." > fragstring_list.txt

for string in `cat fragstring_list.txt`; do

grep -w -A3 "$string" */output.denovo* | head -n3 | tail -n1 >> scorestrings.txt

done

awk -F',' '{for (i=2; i<NF; i++) printf $i ","; print $NF}' scorestrings.txt >>scorestrs.txt
rm scorestrings.txt

cd ..

