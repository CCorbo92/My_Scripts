export system_file=$1
export docking_dir=$2

rm Outcome_AD.txt
rm RMSD_grep_AD_sort.txt
for system in `  cat ${system_file} `
do
echo ${system}

cd ${system}/${docking_dir}

RMSD1=`head -n 1 RMSD_AD_grep.txt  `

sort -n RMSD_AD_grep.txt >> RMSD_grep_AD_sort.txt
RMSD2=`head -n 1 RMSD_grep_AD_sort.txt`

if (( $(echo "$RMSD1 < 2.00" | bc -l) )); then
echo -n "Success " >> ../../Outcome_AD.txt
echo ${system} >> ../../Outcome_AD.txt

elif (( $(echo "$RMSD2 < 2.00" | bc -l) )); then
echo -n "Score Fail " >> ../../Outcome_AD.txt
echo ${system} >> ../../Outcome_AD.txt

else 
echo -n "Sample Fail " >> ../../Outcome_AD.txt
echo ${system} >> ../../Outcome_AD.txt
fi

cd ../../

done
