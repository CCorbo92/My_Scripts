list="countAD.txt"
mkdir Final_Stats_AD
num2=`wc -l Lists_Of_PDB_by_descriptor/List_with_0_7_RB.txt | awk '{print $1}'`
num3=`wc -l Lists_Of_PDB_by_descriptor/List_with_8_15_RB.txt| awk '{print $1}'`
num4=`wc -l Lists_Of_PDB_by_descriptor/List_with_16_up_RB.txt| awk '{print $1}'`
num5=`wc -l Lists_Of_PDB_by_descriptor/List_FDA_PDB.txt| awk '{print $1}'`
num6=`wc -l Lists_Of_PDB_by_descriptor/List_with_ZN.txt| awk '{print $1}'`
for i in `cat ${list}`;
do
#This will print the intersection of file1 and file2
awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Success_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_0_7_RB.txt >> ADSuccess_${i}_List_with_0_7_RB
num1=`wc -l ADSuccess_${i}_List_with_0_7_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num2 " | bc >> ADSuccess_List_with_0_7_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Success_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_8_15_RB.txt >> ADSuccess_${i}_List_with_8_15_RB 
num1=`wc -l ADSuccess_${i}_List_with_8_15_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num3 " | bc >> ADSuccess_List_with_8_15_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Success_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_16_up_RB.txt >> ADSuccess_${i}_List_with_16_up_RB
num1=`wc -l ADSuccess_${i}_List_with_16_up_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num4 " | bc >> ADSuccess_List_with_16_up_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Success_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_FDA_PDB.txt >> ADSuccess_${i}_List_FDA
num1=`wc -l ADSuccess_${i}_List_FDA | awk '{print $1}'`
echo "scale=4 ; $num1 / $num5 " | bc >> ADSuccess_List_FDA_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Success_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_ZN.txt >> ADSuccess_${i}_List_ZN
num1=`wc -l ADSuccess_${i}_List_ZN | awk '{print $1}'`
echo "scale=4 ; $num1 / $num6 " | bc >> ADSuccess_List_ZN_All.txt
#########
awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Score_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_0_7_RB.txt >> ADScore_Fail_${i}_List_with_0_7_RB
num1=`wc -l ADScore_Fail_${i}_List_with_0_7_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num2 " | bc >> ADScore_Fail_List_with_0_7_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Score_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_8_15_RB.txt >> ADScore_Fail_${i}_List_with_8_15_RB
num1=`wc -l ADScore_Fail_${i}_List_with_8_15_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num3 " | bc >> ADScore_Fail_List_with_8_15_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Score_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_16_up_RB.txt >> ADScore_Fail_${i}_List_with_16_up_RB
num1=`wc -l ADScore_Fail_${i}_List_with_16_up_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num4 " | bc >> ADScore_Fail_List_with_16_up_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Score_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_FDA_PDB.txt >> ADScore_Fail_${i}_List_FDA
num1=`wc -l ADScore_Fail_${i}_List_FDA | awk '{print $1}'`
echo "scale=4 ; $num1 / $num5 " | bc >> ADScore_Fail_List_FDA_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Score_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_ZN.txt >> ADScore_Fail_${i}_List_ZN
num1=`wc -l ADScore_Fail_${i}_List_ZN | awk '{print $1}'`
echo "scale=4 ; $num1 / $num6 " | bc >> ADScore_Fail_List_ZN_All.txt
#########
awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Sample_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_0_7_RB.txt >> ADSample_Fail_${i}_List_with_0_7_RB
num1=`wc -l ADSample_Fail_${i}_List_with_0_7_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num2 " | bc >> ADSample_Fail_List_with_0_7_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Sample_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_8_15_RB.txt >> ADSample_Fail_${i}_List_with_8_15_RB
num1=`wc -l ADSample_Fail_${i}_List_with_8_15_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num3 " | bc >> ADSample_Fail_List_with_8_15_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Sample_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_16_up_RB.txt >> ADSample_Fail_${i}_List_with_16_up_RB
num1=`wc -l ADSample_Fail_${i}_List_with_16_up_RB | awk '{print $1}'`
echo "scale=4 ; $num1 / $num4 " | bc >> ADSample_Fail_List_with_16_up_RB_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Sample_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_FDA_PDB.txt >> ADSample_Fail_${i}_List_FDA
num1=`wc -l ADSample_Fail_${i}_List_FDA | awk '{print $1}'`
echo "scale=4 ; $num1 / $num5 " | bc >> ADSample_Fail_List_FDA_All.txt

awk 'FNR==NR{a[$1];next}($1 in a){print}' ./AutoDock_5_Seeds/Sample_Fail_2021_09_09_S${i}_${i}.txt Lists_Of_PDB_by_descriptor/List_with_ZN.txt >> ADSample_Fail_${i}_List_ZN
num1=`wc -l ADSample_Fail_${i}_List_ZN | awk '{print $1}'`
echo "scale=4 ; $num1 / $num6 " | bc >> ADSample_Fail_List_ZN_All.txt

done
mv *All.txt Final_Stats_AD
rm AD*
