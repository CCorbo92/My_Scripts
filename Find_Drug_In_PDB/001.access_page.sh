# This script will access Drugcentral.org for each drug in a list of approved drugs and pull its uniprot code for its 
# validated mechanism of action which can then be used later on to search the PDB for matches of the uniprot code and the drug name
# Written by: Christopher Corbo
# Affiliation: Riz/o Lab, Stony Brook University
# Last Edit Date: 08/2021
# Last Edit by: Christopher Corbo

#Cannot start off with any files 
rm *q=*

cp Approved_Drug_Names.txt >> Approved_Drugs_By_Name.txt
#Converts text to html
sed -i 's/(/%28/g' Approved_Drug_Names.txt
sed -i 's/)/%29/g' Approved_Drug_Names.txt
mv Approved_Drug_Names.txt Approved_Drug_Names_HTML.txt
mv Approved_Drugs_By_Name.txt Approved_Drug_Names.txt

file="Approved_Drug_Names_HTML.txt"

for name in `cat ${file}`; do #Looping on all drugs
txt_name=`echo $name | sed 's/%28/(/g'|  sed 's/%29/)/g'`
wget https://drugcentral.org/?q=${name}
echo -n "drugcentral.org" >> temp_url.txt
grep -m1 "drugcard" index.html\?q\=${txt_name} | cut -c 10-| rev | cut -c 3- | rev >> temp_url.txt
rm index*
while read line;
do
wget ${line}
done < "temp_url.txt"

#Dont need full name of file which changes per drug
mv *q=* temp_page.txt

#This line specifies how many confirmed targets there are which will be used in a for loop for clean consistent output even with varying numbers
num_targ=`grep -B30 "checked" temp_page.txt | grep "uniprot"| tr -d "</>=adefgivhtpunrow:." | sed  's/\"/ /g' | awk '{print $1}' | wc -l`
grep -B30 "checked" temp_page.txt | grep "uniprot"| tr -d "</>=adefgivhtpunrow:." | sed  's/\"/ /g' | awk '{print $1}' >> temp_uniprot.txt
for i in $(seq 1 $num_targ); 
do 
echo -n ${txt_name} "," >> Uniprot_Targets.txt
sed -n "${i}p" temp_uniprot.txt >> Uniprot_Targets.txt
done



rm temp*
done #after looping on all drugs
