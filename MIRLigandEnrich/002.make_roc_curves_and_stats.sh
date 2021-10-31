# This script will read a rank sorted list of actives and decoys and plot ROC Curves and calculate AUC and Enrichment Factor 1% values 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 10/2021
# Last Edit by: Christopher Corbo
system="1UYG"
system_file="${system}_set.txt"
mkdir ${system}
cd ${system}

grep -A1 "MOL" /gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset/zzz.DUDE_Files/${system}/actives_final.mol2 | grep "CHEM" | sort -n | uniq > actives_names.txt

grep -A1 "MOL" /gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset/zzz.DUDE_Files/${system}/decoys_final.mol2 | grep "ZIN" | sort -n | uniq > decoys_names.txt

cd ..

for systemrec in `cat ${system_file}`; do
      
   cd ${system}_${systemrec} 

   rm ${system}_${systemrec}_actives.txt ${system}_${systemrec}_decoys.txt actives_names.txt decoys_names.txt

   echo "Active"
   while read name; do
   # If there are molecules without a single protomer to successfully dock, it will be assigned a 0.000 else it will have its best scoring protomer as score
   if grep -q "${name}" Active_score.txt ; then
      grep ${name} Active_score.txt | sort -n -k2 | head -n1 >> ${system}_${systemrec}_actives.txt
   else
      echo ${name} " 0.00 Active" >> ${system}_${systemrec}_actives.txt
   fi 

   done < "../${system}/actives_names.txt"
   
   echo "Decoy"
   while read name; do

   if grep -q "${name}" Decoy_score.txt ; then
      grep ${name} Decoy_score.txt | sort -n -k2 | head -n1 >> ${system}_${systemrec}_decoys.txt
   else
      echo ${name} " 0.00 Decoy" >> ${system}_${systemrec}_decoys.txt
   fi

   done < "../${system}/decoys_names.txt"
   

   #echo -n ${system} " " >> ../Statistics.txt
   #python ../plot_roc.py ${act_count} ${dec_count} ${system} >> ../Statistics.txt
   #mv ${system}*.png ../plots

   cd ..
done

i=1
while read system_name; do
    declare var${i}=`echo ${system_name}`
    i=$((i+1))
done < "${system}_set.txt"
echo $var1
echo $var2
echo $var3
#I cannot see how to avoid manual input at this step (ie var1 var2 var3)

paste -d " " ${system}_${var1}/${system}_${var1}_actives.txt ${system}_${var2}/${system}_${var2}_actives.txt ${system}_${var3}/${system}_${var3}_actives.txt >> ./${system}/all_rec_name_aligned_active.txt

paste -d " " ${system}_${var1}/${system}_${var1}_decoys.txt ${system}_${var2}/${system}_${var2}_decoys.txt ${system}_${var3}/${system}_${var3}_decoys.txt >> ./${system}/all_rec_name_aligned_decoy.txt
   
