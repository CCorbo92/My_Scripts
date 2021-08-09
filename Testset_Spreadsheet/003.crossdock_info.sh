# This script extracts the mmaker descriptors used in crossdocking family alignment from the chimera.out files and places it into a comma separated file
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 08/2021
# Last Edit by: Christopher Corbo

testset="/gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset"

crossdock_dir="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking"


fam_list="${crossdock_dir}/zzz.sample_lists/family_more_than7.txt"
echo "System,CrossDock_Family,Ca_Atoms_Align,RMSD_to_family_ref" >>CrossDock.csv

for family in `cat ${fam_list}`; do

system_file="${crossdock_dir}/zzz.sample_lists/${family}.txt"

for system in `cat ${system_file}`; do

echo ${system}

matched_atom=`grep "RMSD between" ${crossdock_dir}/Alignment/${family}/${system}/chimera.out| awk '{print $3}'`
RMSD=`grep "RMSD between" ${crossdock_dir}/Alignment/${family}/${system}/chimera.out| awk '{print $8}'`

echo ${system} "," ${family} "," ${matched_atom} "," ${RMSD} >>CrossDock.csv 

done

done
