# This script extracts the DOCK and RDKit descriptors from the mol2 headers and places it into a comma separated file
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 08/2021
# Last Edit by: Christopher Corbo

testset="/gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset"

system_file="clean.systems.all"

echo "Molecular_Weight,DOCK_Rotatable_Bonds,Formal_Charge,HBond_Acceptors,HBond_Donors,Stereocenters,cLogP,TPSA,SA_Score,QED,ESOL,num_of_PAINS,SMILES" >> Descriptors.csv

for system in `cat ${system_file}`; do
echo ${system}
mw=`grep "Molecular_Weight"  ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
rb=`grep "DOCK_Rotatable_Bonds" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
fch=`grep "Formal_Charge" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
hba=`grep "HBond_Acceptors" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
hbd=`grep "HBond_Donors" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
sc=`grep "Stereocenters" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
logp=`grep "cLogP" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
tpsa=`grep "TPSA" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
sa=`grep "SA_Score" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
qed=`grep "QED" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
esol=`grep "ESOL" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
pain=`grep "num_of_PAINS" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `
smile=`grep "SMILES" ${testset}/${system}/${system}.lig.descriptors_scored.mol2 | awk '{print $3}' `

echo $system "," $mw "," $rb "," $fch "," $hba "," $hbd "," $sc "," $logp "," $tpsa "," $sa "," $qed "," $esol "," $pain "," $smile >> Descriptors.csv

done
sort -n Descriptors.csv >> Descriptors_sort.csv
