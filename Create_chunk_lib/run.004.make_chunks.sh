# WRONG VERSION OF THIS SCRIPT GOT COPIED 


# 05-01-17 Jiaye: This script will filter out molecules not assigned with partial charges or with too many rotatable bonds (> 15).
#
Chunk_Size="300000"
Interval_Size="2999"
Library_Size="9399354"
cd ./mol2_decompress/filtered_mol2
mkdir zzz.intervals
i="1"
c="0"
#Bash doesnt allow two variables so set less than value that Count (line below) would produce
#Count=`wc -l Set_A_Lines.txt`
while [ $i -lt 9399354 ]
do
start=`sed -n "${i}p" line_start_stop.txt | awk '{print $1}'` 
echo -n ${start} >> zzz.intervals/Interval${c}.lines.txt 

i=$[$i+${Interval_Size}]

stop=`sed -n "${i}p" line_start_stop.txt | awk '{print $2}'` 
echo "  " ${stop} >> zzz.intervals/Interval${c}.lines.txt
i=$[$i+1]
c=$[$c+1]
echo $c
done

#chunks=`ls -l Set_A/ | grep -v "total" | wc -l`
#replace 11 with num_chunks
#for (( c=0; c<=${chunks}; c++ ))
#do
#   while read start stop;
#   do
#      k=$[$stop+1]
#      sed -n "${start},${stop}p;${k}q" All_Scored.mol2 >> Set_A/chunk${c}.mol2
#   done<"Set_A/chunk${c}.lines.txt"

#   while read start stop;
#   do
#      k=$[$stop+1]
#      sed -n "${start},${stop}p;${k}q" All_Scored.mol2 >> Set_B/chunk${c}.mol2
#   done<"Set_B/chunk${c}.lines.txt"

#   while read start stop;
#   do
#      k=$[$stop+1]
#      sed -n "${start},${stop}p;${k}q" All_Scored.mol2 >> Set_C/chunk${c}.mol2
#   done<"Set_C/chunk${c}.lines.txt"
