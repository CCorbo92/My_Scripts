#!/bin/bash
#The purpose of this software is to automate the pose reproduction with the AutoDock software
#This first arguement will be an input file with the name of all the systems present, an example file format is within /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/clean.systems.all
# Set these LOCALLY IN THE SCRIPT
#This 2nd arguement cd into the main directory where the docking experiment is performed

#module load shared mgltools/1.5.6
#IF YOU HAVE A PROBLEM I BELIEVE THIS WOULD BE THE MOST LIKELY CAUSE
#Seawolve might change these modules and their locations, if this occurs contact seawulf about this issue and ask them where the module was relocated to and then change the script to match the new location
#This module allows you to create the pdbqt files of the receptor and ligand, and generate the grid & docking parameter file
#It's also common to use the AutoDockTools software to visual the molecules manualy, to check your structures
module load mgltools/1.5.6
#This module generates both the grid and performs the actual molecular docking and is essential
module load autodock/4.2.6
#python /gpfs/projects/rizzo/yuchzhou/RCR/DOCK_testset/slaverty_autodock_runs/slaverty_autodock_20_GA_runs_DOCK_grids/AutoGrid_space.py ${systems}
ref_fam=$1
comp_system=$2
dock_dir=$3
list_of_sys="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/${ref_fam}.txt"

for ref_system in `cat ${list_of_sys}`; do
cd ${ref_system}

mkdir ${dock_dir}
cd ${dock_dir}

echo ${ref_system}
#This will be used to create the docking parameter file and to make modification to the docking parameters follow the comments for prepare_gpf4.py
echo "Docking parameter file is being generated with ligand ${comp_system} and receptor ${ref_system}"
/gpfs/software/mgltools/bin/prepare_dpf42.py -l ../${comp_system}.lig.gast.pdbqt -r ../${ref_system}.rec.clean.pdbqt -o ${ref_system}.dock.parameter.dpf
sed -i '2i parameter_file /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/PoseReprod_Autodock/AD4_parameters_with_Na_K.dat # force field default parameter file' ${ref_system}.dock.parameter.dpf
input="${ref_system}.dock.parameter.dpf"
x=0
rm ${ref_system}.docking.dpf
touch ${ref_system}.docking.dpf
     while IFS= read -r line
     do
             #echo "$line"
         map="map"
         if [[ "$line" == *"$map"* ]]; then
             if [[ "$line" == *"fld"* ]]; then
                 echo "It's a field"
                 echo "fld ../${ref_system}.rec.clean.maps.fld # grid_data_file" >> ${ref_system}.docking.dpf
                #sed "${x}d" 1B9V.docking.dpf > 1B9V.docking.dpf
                #sed -i "${x}i fld ../receptor/${lig}/1B9V.rec.clean.fld # grid_data_file"
                #elecmap 1B9V.rec.clean.e.map         # electrostatics map
                #esolvmap 1B9V.rec.clean.d.map       # desolvation map
             elif [[ "$line" == *"elecmap"* ]]; then
                 echo "electrostatic"
                 echo "elecmap ../${ref_system}.rec.clean.e.map # electrostatics map" >> ${ref_system}.docking.dpf
             elif [[ "$line" == *"desolvmap"* ]]; then
                 echo "desolvation"
                 echo "desolvmap ../${ref_system}.rec.clean.d.map # desolvation map" >> ${ref_system}.docking.dpf
             else
                 echo "atom type map"
                    IFS=' ' read -r -a array <<< "$line"
                    echo "map ../${array[1]}  # atom-specific affinity map" >> ${ref_system}.docking.dpf
                    echo ${array[1]}
             fi
         elif [[ "$line" == *"move"* ]]; then
             echo "move ../${comp_system}.lig.gast.pdbqt # small molecule" >> ${ref_system}.docking.dpf
             echo "I DID THIS THING"
         else
             echo "$line" >> ${ref_system}.docking.dpf
         fi
         x=$(( x + 1 ))
     done < "$input"
     sed -i 's/pid time/5       /g' ${ref_system}.docking.dpf
     echo "DOCKING process is occuring"
     grep "move" ${ref_system}.docking.dpf
     autodock4 -p ${ref_system}.docking.dpf -l ${ref_system}.docking.dlg
     echo "DEBUGGING ECHO"
     summarize_results4.py -d /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/21_07_06/zzz.crossdock/${ref_fam}/${comp_system}/${ref_system}/21_07_07_srun 
     #sed -i 's/ga_run 10/ga_run 20/g' ${ref_system}.dock.parameter.dpf
     #Performes the actual autodock experiment
     cd ..
     echo "Directory after first cd"
     pwd
     cd ..
     echo "Directory after second cd"
     pwd
done
