#!/bin/bash

#This script is used to prepare the ligand and receptors to pdbqt
#Receptors in amber charges, the -C option keeps the original charge of the testset system
#The ligand is given gasteiger charges using AutoDock Tools which relies on obabel

#This is the list of systems to be used within the testset
export systems=$1
#This module allows you to create the pdbqt files of the receptor and ligand, and generate the grid & docking parameter file
#It's also common to use the AutoDockTools software to visual the molecules manualy, to check your structures
module load mgltools/1.5.6
#This module generates both the grid and performs the actual molecular docking and is essential
module load autodock/4.2.6
module load anaconda/2
export testset="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.builds"
#This goes through the file containing all of the systems that are being tested given in the first arguement
for SYSTEM in `  cat ${systems} `
#The purpose of this for loop is to both create a directory for each system and then converts them to pdbqt files
do
     cd ${SYSTEM}
     echo ${SYSTEM}
     echo "Preparing Ligands"
     /gpfs/software/mgltools/bin/prepare_ligand4.py -l ${testset}/${SYSTEM}/001.lig-prep/${SYSTEM}.lig.gast.mol2 -o ${SYSTEM}.lig.gast.pdbqt -C
     #/gpfs/software/mgltools/bin/prepare_ligand4.py -l ${testset}/${SYSTEM}/${SYSTEM}.lig.gast.mol2 -U 'nphs' -o ${SYSTEM}.lig.gast.pdbqt 
     #This converts the protein from mol2 to pdbqt in structure and the -C conserves the charges of the ligand
     echo "Prepare Receptor"
     /gpfs/software/mgltools/bin/prepare_receptor4.py -r ${testset}/${SYSTEM}/zzz.autodockfiles/${SYSTEM}.rec.clean.mol2  -o ${SYSTEM}.rec.clean.pdbqt -C 
     cd ..
done
#/gpfs/software/Anaconda3/bin/python /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/slaverty_autodock_runs/slaverty_autodock_Na_manually_11_12_19/Assigning_Ion_Charges_AutoDock4_v3.py ${systems} ${csv}
