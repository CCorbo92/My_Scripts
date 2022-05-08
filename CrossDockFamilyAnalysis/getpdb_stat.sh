# This script pulls the UNIPROT ID from a list of PDB 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 05/2022
# Last Edit by: Christopher Corbo

cd UniProt_Matching

for family in `cat ../family_more_than7.txt`
do

echo ${family}

for sys in `cat ../${family}.txt`
do
#insert pdb link here
wget -q https://www.rcsb.org/structure/${sys}


if grep -q "uniprot" $sys; then
  uniprot=`cat $sys | grep -o -P '.uniprot/.{0,6}' | cut -c 10- | head -n1`
else
  uniprot="n/a"
fi

echo $sys " " $uniprot >> ${family}_stats.txt 


#remove the page
rm ${sys}
done
#Creates a list with only one PDB example of each UniProt ID in a particular family
sort -u  -k2,2 ${family}_stats.txt > ${family}_Uniprot_uniq.txt
done
