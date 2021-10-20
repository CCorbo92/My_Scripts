mkdir splitmol2_temp
mkdir Viable
cd splitmol2_temp
csplit -z --quiet ../21_10_18_AllSelected_edited.mol2 /Name/ '{*}'
grep "viewdock state: V" * | awk -F':' '{print $1}' >> Viable_list.txt
while read line; do
	mv ${line} ../Viable
done<"Viable_list.txt"
