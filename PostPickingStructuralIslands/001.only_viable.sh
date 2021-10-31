rm -r Viable
mkdir splitmol2_temp
mkdir Viable


cd splitmol2_temp

csplit -z --quiet ../21_10_18_OnlyViable_Fullheader.mol2 /Name/ '{*}'

grep "viewdock state: V" * | awk -F':' '{print $1}' > Viable_list.txt

while read line; do
	mv ${line} ../Viable
done<"Viable_list.txt"

cd ../Viable
grep "Hungarian" * |awk '{print $3 " " $1}' | sort -n |  awk -F':' '{print $1}'| awk '{print $2}'  > Sorted_HMS_Viable.txt

mkdir Rank_Sorted_Mol2

cd Rank_Sorted_Mol2
while read line; do

cat ../$line >> Viable_HMS_sorted.mol2

done<"../Sorted_HMS_Viable.txt"

mkdir Viable_Sorted
cd Viable_Sorted
csplit -z --quiet  ../Viable_HMS_sorted.mol2 /Name/ '{*}'

ls -l | tr -d 'x'| awk '{print $9}' | sort -n > ../../../viable_sorted.txt
