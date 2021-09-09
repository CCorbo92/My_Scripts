
# 05-01-17 Jiaye: This script will filter out molecules not assigned with partial charges or with too many rotatable bonds (> 15).
#

cd ./mol2_decompress/filtered_mol2

Count=`ls -l | grep "scored" | wc -l`
echo $Count
i="0"
#Cat all the separate mol2 into a single one (Will not be okay for very large libraries)
while [ $i -lt ${Count} ]
do
cat ${i}.output_scored.mol2 >> All_Scored.mol2
i=$[$i+1]
done
#Find the start line for each individual molecule in a mol2
grep -n "Name:" All_Scored.mol2 | awk -F':' '{print $1}' >> line_start.txt
grep -n "TEMP" All_Scored.mol2 | awk -F':' '{print $1}' >> line_stop.txt
paste -d ' ' line_start.txt line_stop.txt >> line_start_stop.txt

rm line_start.txt line_stop.txt


