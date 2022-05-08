# This script 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: --/----
# Last Edit by: Christopher Corbo
for family in `cat /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/family_more_than7.txt`
do
 echo $family
 cd $family
 rm PercentSimilarity.txt *.fasta 
 for file in `cat returned_list.txt`;
 do
   cat $file >> ${family}.fasta
   python ../align.py $family
 done
 python ../calculate_similarity.py $family >> PercentSimilarity.txt
 cd ..
done
