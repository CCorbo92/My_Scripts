#!/bin/tcsh -fe

#
# This script will postprocess the output files from virtual screening on cluster. Briefly, it
# concatenates everything together, generates a ranked csv file consisting of ZINC ids and other
# descriptors, calculates some additional descriptors in MOE, clusters by fingerprint, and writes
# some final output files ranked by different scoring functions. It is these files that should be
# visually inspected for purchasing.
#
# Prior to running this script, set the max number of molecules that should be clustered. If the
# number is set to around 100,000 then you can expect the calculation to take about a week. (?)
#
module load anaconda/3

set max_num = "${MAX_NUM_MOL}"


### Set some paths
set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"


set wcl   = 12:00:00
set nodes = 1
set ppn   = 24
set queue = "long-24core"
@ numprocs = (${nodes} * ${ppn})

# If doing best first clustering, this should be set to "_bestfit", else ""
set type = "_bestfit"
set sim  = "0.95"
#or
# set type = ""
# set sim  = "0.75"

### Make a directory for compiling all of the docked results for a given vendor. If the top
### directory already exists, don't remove other vendor results.
if (! -e ${rootdir}/${system}/010.rdkit-and-cluster/) then
	mkdir -p ${rootdir}/${system}/010.rdkit-and-cluster/
endif

if (! -e ${rootdir}/${system}/010.rdkit-and-cluster/${vendor}) then
	mkdir -p ${rootdir}/${system}/010.rdkit-and-cluster/${vendor}
endif
# This is how the naming of the scores was accomplished in the last step so it simplifies it to keep with this nomenclature

foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)
  echo ${score}

  rm -rf ${rootdir}/${system}/010.rdkit-and-cluster/${vendor}/${score}_rank
  mkdir -p ${rootdir}/${system}/010.rdkit-and-cluster/${vendor}/${score}_rank
  cd ${rootdir}/${system}/010.rdkit-and-cluster/${vendor}/${score}_rank

### Count the number of chunks
set num_chunks = `ls -l ${rootdir}/${system}/005.zinclibs/${vendor}/ | grep chunk | wc -l`
echo "num_chunks = ${num_chunks}"
echo "Submitting to queue..."




### Copy the relevant files over
cp ${rootdir}/${system}/009.rank_per_score/${vendor}/${score}_top_${max_num}.mol2 ${system}.${vendor}.${score}_rank.total.mol2
cp ${rootdir}/${system}/009.rank_per_score/${vendor}/RDKit_FP.output_Continuous_Score:_footprint_scored.txt ${system}.${vendor}.${score}_rank.total_fp.txt

# Need footprint text fileCPC
### 1. Check for any duplicate ZINC names, and only keep the molecule with the best Footprint_Similarity_Score.
echo "Checking for duplicate ZINC names..."
# THIS CAN BE DONE ON THE CSV IN 009 by printing out the duplicate zinc ids (This txt file does not actually get used but can be checked for sanity reasons
awk -F',' '{print $1}' |  sort -n  ${rootdir}/${system}/009.rank_per_score/${vendor}/${score}_sorted_${max_num}.csv| uniq -d >> List_duplicates.txt

### 2. Perform a best first clustering on the MACCS fingerprints
echo "Running best first clustering"
set file = "${rootdir}/${system}/009.rank_per_score/${vendor}/ZINC_MACCS_${score}.txt"

python ${scriptdir}/best_first_maccs_clustering.py ${file} >> cluster_py.out
# More straightforward to process the output of the python script as post processing but could be made faster if the output could 
# be manipulated more readily within the python script
set i = 1
#max_num + 1
while ($i < 100001)
set line="`sed -n "${i}p" ${file}`"
set ZName=`echo ${line} | awk -F',' '{print $1}'`
set CName=`grep -B1 ${ZName}  cluster_py.out | head -n1 | awk '{print $3}'`
set Size=`grep -A1 ${ZName}  cluster_py.out | tail -n1`
echo  ${ZName} "," ${CName} "," ${Size} >> Cluster_Info_${score}.txt
   @ i++
end

#
### 6. Combine the csvs containing MOE descriptors and DOCK scores, and print out csv files of just
###    the clusterheads, and of everything grouped by cluster into families 



echo "Combining DOCK / RDKIT descriptors and Cluster Info..."
paste -d ',' Cluster_Info_${score}.txt  ${rootdir}/${system}/009.rank_per_score/${vendor}/${score}_sorted_100000.csv >> Descriptors_${score}_${max_num}.csv
# This will sort based on the cluster number and make unique on the clusters so that only the clusterheads will be in this list
# It was already sorted based on the primary score such that the clusterhead will be the highest on the primary score
sort -u -n -t, -k2 Descriptors_${score}_${max_num}.csv >> Clusterheads_${score}_only.csv

# This step takes the csv of 100000 molecules and removes the lower rank duplicate ZINC ID
sort -u -k1,1 Descriptors_${score}_${max_num}.csv >> Descriptors_${score}_${max_num}_Unique.csv

set file = "Clusterheads_${score}_only.csv"
#Printing out the top 1000 clusterheads for each secondary score
# This line will sort the csv of clusterheads for each secondary score and take the top 1000 ranked
python ${scriptdir}/sorting_scores.py ${file} >> clusterheads_${score}_primary.out
foreach second_score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: )
grep -A1000 ${second_score} clusterheads_${score}_primary.out | tail -n1000 >> clusterheads_${second_score}_secondary.csv

echo "Secondary "  ${second_score}
# This step will determine the cluster size for each of the top 1000 clusterheads per secondary score
awk -F, '{print $3}' clusterheads_${second_score}_secondary.csv >> clusterheads_${second_score}_secondary_cluster_sizes.txt
#Remove leading spaces that are a carry over
sed -i "s/  //g"  clusterheads_${second_score}_secondary_cluster_sizes.txt
#DO a for loop on i and exit once sum is greater than 1000
#This step determines how many clusters it takes to reach 1000 for the family mol2s
set temp_sum=0
set chunk_counter=1
while ($temp_sum < 1000)
   set temp_sum=`awk -v a=$chunk_counter 'NR==1,NR==a {sum+=$1;} END{print sum;}' clusterheads_${second_score}_secondary_cluster_sizes.txt`
   #echo "temp sum " $temp_sum
   @ chunk_counter++
end
#CAN DELETE THIS LINE BELOW LATER ON
echo  ${second_score} "," ${chunk_counter} >> Families_Counting.txt
#Once you know how many clusters to include for the families list then print out the cluster name for each of these clusters
awk -F, '{print $2}' clusterheads_${second_score}_secondary.csv | head -n ${chunk_counter} >> Family_Clusters_${second_score}_secondary.txt
sed -i "s/  //g" Family_Clusters_${second_score}_secondary.txt
#grep each element from clusterheads_${second_score}_secondary_cluster_list.txt  
#You can use this line once you have the clusters # that make it into the family mol2
#Search each of the cluster names in the CSV of 100,000 molecules made unique and pulls out the line for each molecule in the cluster to include
set length=`wc -l Family_Clusters_${second_score}_secondary.txt | awk '{print $1}'`
set i = 1
while ($i < $length)
   set line=`sed -n ${i}p Family_Clusters_${second_score}_secondary.txt`
   awk -v c=$line -F, '{if ($2 == c ) print $0;}' Descriptors_${score}_${max_num}_Unique.csv >> families_${second_score}_secondary.csv
   @ i++
end


#DO FOR ELEMENT IN Family_Clusters_${second_score}_secondary.txt REPLACE 30226 WITH ELEMENT IN LINE BELOW
#awk -F, '{if ($2 == 30226 ) print $0;}' clusterheads_Continuous_Score\:_secondary.out

end #End secondary scoring loop

echo "Job finished"



### Submit the script
#qsub ${system}.${vendor}.${score}.postprocess.qsub.csh

end
exit


