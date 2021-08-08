# This script will access Drugcentral.org for each drug in a list of approved drugs and pull its uniprot code for its 
# validated mechanism of action which can then be used later on to search the PDB for matches of the uniprot code and the drug name
# Written by: Christopher Corbo
# Affiliation: Riz/o Lab, Stony Brook University
# Last Edit Date: 08/2021
# Last Edit by: Christopher Corbo
cp Approved_Drug_Names.txt >> Approved_Drugs_By_Name.txt
sed -i 's/(/%28/g' Approved_Drug_Names.txt
sed -i 's/)/%29/g' Approved_Drug_Names.txt
mv Approved_Drug_Names.txt Approved_Drug_Names_HTML.txt

#file="Approved_Drug_Names_HTML.txt"
file="test.txt"

for name in `cat ${file}`; do
txt_name=`echo $name | sed 's/%28/(/g'|  sed 's/%29/)/g'`
wget https://drugcentral.org/?q=${name}
echo -n "drugcentral.org" >> temp_url.txt
grep -m1 "drugcard" index.html\?q\=${txt_name} | cut -c 10-| rev | cut -c 3- | rev >> temp_url.txt
while read line;
do
wget ${line}
done < "temp_url.txt"
cat temp_url.txt | cut -c 26- >> temp_url2.txt
rm temp*
rm *q=* 
done
