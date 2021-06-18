#!/bin/bash

#This script is used to prepare the ligand and receptors to pdbqt
#Receptors in amber charges, the -C option keeps the original charge of the testset system
#The ligand is given gasteiger charges using AutoDock Tools which relies on obabel

#This is the list of systems to be used within the testset
export systems=$1
#This is the directory that prepares each of these systems
export new_dir=$2
#This module allows you to create the pdbqt files of the receptor and ligand, and generate the grid & docking parameter file
#It's also common to use the AutoDockTools software to visual the molecules manualy, to check your structures
module load mgltools/1.5.6
#This module generates both the grid and performs the actual molecular docking and is essential
module load autodock/4.2.6
export testset="/gpfs/scratch/ccorbo/Benchmarking_and_Validation/CrossDocking/zzz.builds/zzz.master"
cd ${new_dir}
#This goes through the file containing all of the systems that are being tested given in the first arguement
for SYSTEM in `  cat ${systems} `
#The purpose of this for loop is to both create a directory for each system and then converts them to pdbqt files
do
     cd ${SYSTEM}
     echo ${SYSTEM}
     #This creates a working directory for each system
     #This prepares the ligand converting the structures from mol2 files to pdbqt which autodock uses
     #The -U '' in option in prepare_ligand4.py is used to clean up ligand and the '' means to make little to no changes in the structure
     #Default is to add gasteiger charges to the molecule
     #IMPORTANT NOTE, SHOULD READ IF EXPERIENCING ERRORS, REGARDING LOCATION OF ORIGINAL  mol2 files
     #This code takes from yuchzhou directory of the systems,so if this system isn't updated with the necessary systems and files then this code will likely not work
     #Solution to this issue, I'd recommend to change the location that the original mol2 is being copied from to solve this issue
     echo "Preparing Ligands"
     /gpfs/software/mgltools/bin/prepare_ligand4.py -l ${testset}/${SYSTEM}.lig.moe.gast.mol2 -o ${SYSTEM}.lig.gast.pdbqt -C
     #/gpfs/software/mgltools/bin/prepare_ligand4.py -l ${testset}/${SYSTEM}/${SYSTEM}.lig.gast.mol2 -U 'nphs' -o ${SYSTEM}.lig.gast.pdbqt 
     #This converts the protein from mol2 to pdbqt in structure and the -C conserves the charges of the ligand
     echo "Prepare Receptor"
     /gpfs/software/mgltools/bin/prepare_receptor4.py -r ${testset}/${SYSTEM}.rec.foramber.mol2  -o ${SYSTEM}.rec.clean.pdbqt -C 
     cd ..
done
#/gpfs/software/Anaconda3/bin/python /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/slaverty_autodock_runs/slaverty_autodock_Na_manually_11_12_19/Assigning_Ion_Charges_AutoDock4_v3.py ${systems} ${csv}
