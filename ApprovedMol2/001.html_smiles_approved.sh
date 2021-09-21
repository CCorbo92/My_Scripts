# This script gets a csv of smiles and names of approved drugs from drugcentral.org and converts 
# the smiles to a searchable html format 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: --/----
# Last Edit by: Christopher Corbo
wget https://unmtid-shinyapps.net/download/structures.smiles.tsv

i="1"
while read -r smile col2 col3 col4 name col6;
do
echo $smile  " " $name
echo ${i}
echo -n ${name} >> smile_html.txt
echo -n " " >> smile_html.txt
echo -n ${smile} >> smile_html.txt
echo -n " " >> smile_html.txt
echo $smile |sed 's/(/%28/g' | sed 's/=/%3D/g'| sed 's/)/%29/g'| sed 's/\[/%5B/g'| sed 's/]/%5D/g'| sed 's/@/%40/g'| sed 's/+/%2B/g'| sed 's/#/%23/g' >> smile_html.txt




i=$(($i+1))
done <"structures.smiles.tsv"

grep -v "SMILE" smile_html.txt >> tmp.txt
rm smile_html.txt
mv tmp.txt smile_html.txt
