#!/bin/bash
#The purpose of this software is to automate the pose reproduction with the AutoDock software
#This script relies on primarily 
export WORKDIR="/gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/zzz.subset.lists"
export testset="/gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset"
#Majority of the results of the AutoDock testset was performed within this directory
export datafolder="/gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/slaverty_autodock_runs"
#This first arguement will be an input file with the name of all the systems present, an example file format is within /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/clean.systems.all
#This inputs the names of the systems that you are working with in your program, an example file is /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/clean.systems.all
export SYSTEM=$1
#This 2nd arguement cd into the main directory where the docking experiment is performed
export new_dir=$2
#This 3rd arguement is the dock directory, this allows you to reuse the receptors, ligands, and grids for all docking 
export dock_dir=$3
#module load shared mgltools/1.5.6
#IF YOU HAVE A PROBLEM I BELIEVE THIS WOULD BE THE MOST LIKELY CAUSE
#Seawolve might change these modules and their locations, if this occurs contact seawulf about this issue and ask them where the module was relocated to and then change the script to match the new location
#This module allows you to create the pdbqt files of the receptor and ligand, and generate the grid & docking parameter file
#It's also common to use the AutoDockTools software to visual the molecules manualy, to check your structures
module load mgltools/1.5.6
#This module generates both the grid and performs the actual molecular docking and is essential
module load autodock/4.2.6
#python /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/slaverty_autodock_runs/slaverty_autodock_20_GA_runs_DOCK_grids/AutoGrid_space.py ${systems}
cd ${new_dir}
     echo ${SYSTEM}
     cd ${SYSTEM}
     mkdir ${dock_dir}
     cd ${dock_dir}
     #This will be used to create the docking parameter file and to make modification to the docking parameters follow the comments for prepare_gpf4.py
     echo "Docking parameter file is being generated"
     /gpfs/software/mgltools/bin/prepare_dpf42.py -l ../${SYSTEM}.lig.gast.pdbqt -r ../${SYSTEM}.rec.clean.pdbqt -o ${SYSTEM}.dock.parameter.dpf
     sed -i '2i parameter_file /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/AD4_parameters_with_Na_K.dat # force field default parameter file' ${SYSTEM}.dock.parameter.dpf
     input="${SYSTEM}.dock.parameter.dpf"
     x=0
     rm ${SYSTEM}.docking.dpf
     touch ${SYSTEM}.docking.dpf
     while IFS= read -r line
     do
             #echo "$line"
         map="map"
         if [[ "$line" == *"$map"* ]]; then
             if [[ "$line" == *"fld"* ]]; then
                 echo "It's a field"
                 echo "fld ../${SYSTEM}.rec.clean.maps.fld # grid_data_file" >> ${SYSTEM}.docking.dpf
                #sed "${x}d" 1B9V.docking.dpf > 1B9V.docking.dpf
                #sed -i "${x}i fld ../receptor/${lig}/1B9V.rec.clean.fld # grid_data_file"
                #elecmap 1B9V.rec.clean.e.map         # electrostatics map
                #esolvmap 1B9V.rec.clean.d.map       # desolvation map
             elif [[ "$line" == *"elecmap"* ]]; then
                 echo "electrostatic"
                 echo "elecmap ../${SYSTEM}.rec.clean.e.map # electrostatics map" >> ${SYSTEM}.docking.dpf
             elif [[ "$line" == *"desolvmap"* ]]; then
                 echo "desolvation"
                 echo "desolvmap ../${SYSTEM}.rec.clean.d.map # desolvation map" >> ${SYSTEM}.docking.dpf
             else
                 echo "atom type map"
                    IFS=' ' read -r -a array <<< "$line"
                    echo "map ../${array[1]}  # atom-specific affinity map" >> ${SYSTEM}.docking.dpf
                    echo ${array[1]}
             fi
         elif [[ "$line" == *"move"* ]]; then
             echo "move ../${SYSTEM}.lig.gast.pdbqt # small molecule" >> ${SYSTEM}.docking.dpf
         else
             echo "$line" >> ${SYSTEM}.docking.dpf
         fi
         x=$(( x + 1 ))
     done < "$input"
     sed -i 's/pid time/5       /g' ${SYSTEM}.docking.dpf
     echo "DOCKING process is occuring"
     autodock4 -p ${SYSTEM}.docking.dpf -l ${SYSTEM}.docking.dlg
     summarize_results4.py -d ../${dock_dir}
     #sed -i 's/ga_run 10/ga_run 20/g' ${SYSTEM}.dock.parameter.dpf
     #Performes the actual autodock experiment
     cd ..
     cd ..

