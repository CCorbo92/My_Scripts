#!/bin/bash

export systems=$1
export new_dir=$2
#This module allows you to create the pdbqt files of the receptor and ligand, and generate the grid & docking parameter file
#It's also common to use the AutoDockTools software to visual the molecules manualy, to check your structures
module load mgltools/1.5.6
module load autodock/4.2.6
echo "Grid parameter file is being generated"
cd ${new_dir}
for SYSTEM in `  cat ${systems} `
do
     cd ${SYSTEM}
     echo ${SYSTEM}
     #This step is creates the grid parameters used to create the grid.
     #The -y is critically important for this project because it tells the software to generate centered around the ligand so it knows where to dock
     #IF YOU WANT TO CHANGE THE GRID PARAMETERS, Follow instructions with these two main methods
     #The npts determines grid size npts='x,y,z', default is npts='40,40,40'
     #The first and preferred way is to input the parameters into prepare_gpf4.py which would require you to read the instructions about this
     #The second is to use the second command to change it which is shown in thee sed command currently commented out below
     #SIDENOTE: The reason the mod.${SYSTEM}.rec.clean.pdbqt is used instead of the previous ${SYSTEM}.rec.clean.pdbqt is because the prepare_receptor4.py was unable to give charges to ions so the Assigning_Ion_Charges_AutoDock4_v2.py program is used to generate these receptors with the proper charges on thm to resolve this issue and as a by product the new mod.${SYSTEM}.rec.clean.pdbqt is used instead
     echo "Grid parameter file is being generated"
     /gpfs/software/mgltools/bin/prepare_gpf4.py -l ${SYSTEM}.lig.gast.pdbqt -r ${SYSTEM}.rec.clean.pdbqt -y #-p npts='60,60,60' -p spacing="0.300" -y
     sed -i '2i parameter_file /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/AD4_parameters_with_Na_K.dat # force field default parameter file' ${SYSTEM}.rec.clean.gpf
     cd ..
done
#module load anaconda/3
#/gpfs/software/Anaconda3/bin/python /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/slaverty_autodock_runs/slaverty_autodock_Na_manually_11_12_19/Assigning_Ion_Charges_AutoDock4_v3.py ${systems} ${csv}
#module unload anaconda/3
for SYSTEM in `  cat ${systems} `
do
     cd ${SYSTEM}
     echo "Grid is being generated"
     autogrid4 -p ${SYSTEM}.rec.clean.gpf -l ${SYSTEM}.autogrid.glg
     cd ..
done

