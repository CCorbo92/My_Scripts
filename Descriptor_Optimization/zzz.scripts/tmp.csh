foreach score (Continuous_Score Pharmacophore_Score Hungarian_Matching_Similarity_Score Property_Volume_Score Footprint_Similarity_Score Tanimoto_Score)
  echo  ${score}_actives.txt
  echo ${score}_decoys.txt

#paste -d ${score}_actives.txt MACCSfinal_\${score}.txt -d "," >
end
