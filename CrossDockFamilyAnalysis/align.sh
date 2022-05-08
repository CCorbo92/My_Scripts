# This script reads in a file previously made *_Uniprot_uniq.txt and then fetches the FASTA of each PDB ID and then aligns them in pairwise alignments. The percent similarity is calculated.
# Need to fix calculate_similarity.py which assumes FASTA heading only has one set of parentheses
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 05/2022
# Last Edit by: Christopher Corbo
for family in `cat /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/family_more_than7.txt`
do
 echo $family
 rm -r $family
 mkdir $family
 cd $family
 while read PDB uniprot; do
 
   wget -q https://www.rcsb.org/fasta/entry/${PDB}/display
 
 done<"/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/UniProt_Matching/${family}_Uniprot_uniq.txt"
 
 tmp=`ls -l | awk '{print $9}'` 
 echo $tmp >>returned_list.txt
 for file in `cat returned_list.txt`;
 do
   cat $file >> ${family}.fasta
   python ../align.py $family
 done
 python ../calculate_similarity.py $family >> PercentSimilarity.txt
 cd ..
done
