#!/bin/bash

#This module allows you to create the pdbqt files of the receptor and ligand, and generate the grid & docking parameter file
#It's also common to use the AutoDockTools software to visual the molecules manualy, to check your structures
module load anaconda/3
module load mgltools/1.5.6
module load autodock/4.2.6

ref_fam=$1
comp_system=$2

list_of_sys="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/${ref_fam}.txt"

for ref_system in `cat ${list_of_sys}`; do
mkdir ${ref_system}
cd ${ref_system}

echo "Grid parameter file is being generated for Pair ${comp_system} and ${ref_system}"
     #This step is creates the grid parameters used to create the grid.
     #The -y is critically important for this project because it tells the software to generate centered around the ligand so it knows where to dock
     #IF YOU WANT TO CHANGE THE GRID PARAMETERS, Follow instructions with these two main methods
     #The npts determines grid size npts='x,y,z', default is npts='40,40,40'
     #The first and preferred way is to input the parameters into prepare_gpf4.py which would require you to read the instructions about this
     #The second is to use the second command to change it which is shown in thee sed command currently commented out below
     #SIDENOTE: The reason the mod.${SYSTEM}.rec.clean.pdbqt is used instead of the previous ${SYSTEM}.rec.clean.pdbqt is because the prepare_receptor4.py was unable to give charges to ions so the Assigning_Ion_Charges_AutoDock4_v2.py program is used to generate these receptors with the proper charges on thm to resolve this issue and as a by product the new mod.${SYSTEM}.rec.clean.pdbqt is used instead
     echo "Grid parameter file is being generated"
     cp /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/21_07_06/${comp_system}/${comp_system}.lig.gast.pdbqt ./
     cp /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/21_07_06/${ref_system}/${ref_system}.rec.clean.pdbqt ./
     /gpfs/software/mgltools/bin/prepare_gpf4.py -l ${comp_system}.lig.gast.pdbqt -r ${ref_system}.rec.clean.pdbqt -y #-p npts='60,60,60' -p spacing="0.300" -y
     sed -i '2i parameter_file /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/PoseReprod_Autodock/AD4_parameters_with_Na_K.dat # force field default parameter file' ${ref_system}.rec.clean.gpf
#module load anaconda/3
#/gpfs/software/Anaconda3/bin/python /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/slaverty_autodock_runs/slaverty_autodock_Na_manually_11_12_19/Assigning_Ion_Charges_AutoDock4_v3.py ${systems} ${csv}
#module unload anaconda/3
     echo "Grid is being generated"
     autogrid4 -p ${ref_system}.rec.clean.gpf -l ${ref_system}.autogrid.glg
     cd ..
done

