file="Approved_Drug_Names.txt"

for drug in `cat ${file}`; do

smiles=`grep -w "${name}" smiles.tsv  | awk '{print $1}'`

echo $smiles

done
