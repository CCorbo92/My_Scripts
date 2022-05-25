#!/bin/tcsh
#SBATCH --partition=rn-long-40core
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=rscor
#SBATCH --output=rscor.out

set system_in = "1B9V"

cd /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/${system_in}

foreach score (Continuous_Score Pharmacophore_Score Hungarian_Matching_Similarity_Score Property_Volume_Score Footprint_Similarity_Score Tanimoto_Score)
  grep $score ${system_in}_Actives_rescore_ranked.mol2 | awk '{print $3}'  > ${score}_actives.txt
  grep $score ${system_in}_Decoys_rescore_ranked.mol2 | awk '{print $3}' > ${score}_decoys.txt

#paste -d ${score}_actives.txt MACCSfinal_\${score}.txt -d "," >
end
cp Continuous_Score_actives.txt Continuous_Score_actives_W.txt
cp Continuous_Score_decoys.txt Continuous_Score_decoys_W.txt

awk  '{ printf("%.3f\n", $0 *(4) ) }' Pharmacophore_Score_actives.txt > Pharmacophore_Score_actives_W.txt
awk  '{ printf("%.3f\n", $0 *(4) ) }' Pharmacophore_Score_decoys.txt > Pharmacophore_Score_decoys_W.txt

awk  '{ printf("%.3f\n", $0 *(3) ) }' Hungarian_Matching_Similarity_Score_actives.txt > Hungarian_Matching_Similarity_Score_actives_W.txt
awk  '{ printf("%.3f\n", $0 *(3) ) }' Hungarian_Matching_Similarity_Score_decoys.txt > Hungarian_Matching_Similarity_Score_decoys_W.txt

awk  '{ printf("%.3f\n", $0 *(-80) ) }' Property_Volume_Score_actives.txt > Property_Volume_Score_actives_W.txt
awk  '{ printf("%.3f\n", $0 *(-80) ) }' Property_Volume_Score_decoys.txt > Property_Volume_Score_decoys_W.txt

awk  '{ printf("%.3f\n", $0 *(2) ) }' Footprint_Similarity_Score_actives.txt > Footprint_Similarity_Score_actives_W.txt
awk  '{ printf("%.3f\n", $0 *(2) ) }' Footprint_Similarity_Score_decoys.txt > Footprint_Similarity_Score_decoys_W.txt

awk  '{ printf("%.3f\n", $0 *(-80) ) }' Tanimoto_Score_actives.txt > Tanimoto_Score_actives_W.txt
awk  '{ printf("%.3f\n", $0 *(-80) ) }' Tanimoto_Score_decoys.txt > Tanimoto_Score_decoys_W.txt
