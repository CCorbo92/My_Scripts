# WRONG VERSION OF THIS SCRIPT GOT COPIED 


# 05-01-17 Jiaye: This script will filter out molecules not assigned with partial charges or with too many rotatable bonds (> 15).
#
Chunk_Size="300000"
Interval_Size="2999"
Library_Size=" 9,399,000"
#Library Size / Chunk Size (Celing Value)
Num_Chunks="31"

cd ./mol2_decompress/filtered_mol2/zzz.intervals
i="0"
c="0"
#Bash doesnt allow two variables so set less than value that Count (line below) would produce
Count=`ls -l | wc -l `
#Count="3163"
cd ..
mkdir zzz.chunks_final

while [ "$c" -lt "$Count" ]; do
echo $c
if [ "$i" -eq "$Num_Chunks" ];
then
i="0"
fi

start=`awk '{print $1}' ./zzz.intervals/Interval${c}.lines.txt `
stop=`awk '{print $2}'  ./zzz.intervals/Interval${c}.lines.txt`
k=$[$stop+1]
sed -n "${start},${stop}p;${k}q" All_Scored.mol2 >> zzz.chunks_final/chunk${i}.mol2

i=$[$i+1]
c=$[$c+1]

done

#Do leftovers
#c=${Count}
#Count=`ls -l ./zzz.intervals | wc -l `

#while [ "$c" -lt "$Count" ]; do
#echo $c

#start=`awk '{print $1}' ./zzz.intervals/Interval${c}.lines.txt `
#stop=`awk '{print $2}'  ./zzz.intervals/Interval${c}.lines.txt`
#k=$[$stop+1]
#sed -n "${start},${stop}p;${k}q" All_Scored.mol2 >> zzz.chunks_final/chunk_leftover.mol2

#c=$[$c+1]

#done

