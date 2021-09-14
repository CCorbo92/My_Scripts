i="1"
while read -r name smile html;
do
echo ${i} 

wget http://zinc.docking.org/substances/search/?q=${html} 

mv index* out.txt
ZINC="$(grep '<a href="/substances/ZINC' out.txt | head -n1| cut -c 30- | rev | cut -c 3-| rev)"
echo ${ZINC}
wget http://zinc.docking.org/substances/${ZINC}/

#link="$(grep '<a href="/substances.mol2?q=' $$$$$ | head -n1 | cut -c 43- | rev | cut -c 3- | rev)"
#wget ${link}
link="$(grep '<li><a href="http://files.docking.org/protomers' index.html | grep "Mol2" | head -n1 | cut -c 38- | rev | cut -c 16- | rev)"
wget ${link}
mv *.mol2.gz ./ZINC/${name}.mol2.gz
rm out.txt
rm index.html
rm *.mol2
i=$(($i+1))

done <"smile_html.txt"

