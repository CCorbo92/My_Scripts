system_file="./clean.systems.all"
testset="/gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset"

mkdir temp

#These correspond to each column in the merged spreadsheet
COL1="System"
COL2="MW"
COL3="DOCK_Rotatable_Bonds"
COL4="Formal_Charge"
COL5="HBondAcceptors"
COL6="HBondDonors"
COL7="Stereocenters"
COL8="cLogP"
COL9="TPSA"
COL10="SA"
COL11="QED"
COL12="ESOL"
COL13="num_of_PAINS"
COL14="SMILES"
COL15="System"
COL16="Drug_Name"
COL17="Year_FDA_Approved"
COL18="System"
COL19="Protein_Type"
COL20="Protein_Subtype"
COL21="System"
COL22="Crossdock_Family"
COL23="Ca_Atoms_Align"
COL24="RMSD_to_family_ref"

for system in `cat ${system_file}`; do  
echo ${system}
name=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $1}'`
mw=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $2}'`
rb=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $3}'`
charge=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $4}'`
hba=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $5}'`
hbd=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $6}'`
stereo=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $7}'`
logp=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $8}'`
tpsa=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $9}'`
sa=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $10}'`
qed=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $11}'`
esol=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $12}'`
pain=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $13}'`
smile=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $14}'`
drugname=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $16}'`
year=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $17}'`
protein=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $19}'`
subtype=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $20}'`
cdfam=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $22}'`
camatch=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $23}'`
rmsd=`grep ${system} Full_spreadsheet.csv | awk -F',' '{print $24}'`

printf "########## %-20s  %-20s" $COL1 $name >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL2 $mw >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL3 $rb >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL4 $charge >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL5 $hba >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL6 $hbd >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL7 $stereo >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL8 $logp >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL9 $tpsa >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL10 $sa >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL11 $qed >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL12 $esol >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL13 $pain >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL14 $smile >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL16 $drugname >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL17 $year >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL19 $protein >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL20 $subtype >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL22 $cdfam >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL23 $camatch >> temp/${system}.header
echo " " >> temp/${system}.header
printf "########## %-20s  %-20s" $COL24 $rmsd >> temp/${system}.header
echo " " >> temp/${system}.header

cat temp/${system}.header >> All_testset.mol2
echo " " >> All_testset.mol2
cat  ${testset}/${system}/${system}.lig.am1bcc.mol2 >> All_testset.mol2


done

rm -r temp
