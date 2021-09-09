
#
cd ./mol2_decompress/filtered_mol2/zzz.chunks_final
mv ../All_Names.txt ../All_Names_old.txt

c="0"
#Bash doesnt allow two variables so set less than value that Count (line below) would produce
#Count=`wc -l Set_A_Lines.txt`

#replace 11 with num_chunks
for (( c=0; c<=39; c++ ))
do
echo ${c}
      grep "Name" chunk${c}.mol2 | awk '{print $3}' >> ../All_Names.txt
done
cd ..
sort -n All_Names.txt | uniq | wc -l 
