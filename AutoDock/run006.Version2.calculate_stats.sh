export system_file=$1
export docking_dir=$2

rm Outcome.txt

for system in `  cat ${system_file} `
do
echo ${system}

cd ${system}/${docking_dir}

Index=`cat index_lowest_score.txt`

RMSD1=`head -n ${Index} RMSD_grep.txt | tail -n1 `

sort -n RMSD_grep.txt >> RMSD_grep_sort.txt
RMSD2=`head -n 1 RMSD_grep_sort.txt`

if (( $(echo "$RMSD1 < 2.00" | bc -l) )); then
echo -n "Success " >> ../../Outcome.txt
echo ${system} >> ../../Outcome.txt

elif (( $(echo "$RMSD2 < 2.00" | bc -l) )); then
echo -n "Score Fail " >> ../../Outcome.txt
echo ${system} >> ../../Outcome.txt

else 
echo -n "Sample Fail " >> ../../Outcome.txt
echo ${system} >> ../../Outcome.txt
fi

cd ../../

done
