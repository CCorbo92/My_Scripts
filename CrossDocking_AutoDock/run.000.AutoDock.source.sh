#!/bin/bash

#Arguement 1 is the list of PDB codes
export systems=$1
#Argument 2 is what to name the new directory to run the docking
export new_dir=$2
mkdir ${new_dir}
#This goes through the file containing all of the systems that are being tested given in the first arguement
for SYSTEM in `  cat ${systems} `
#The purpose of this for loop is to both create a directory for each system and then converts them to pdbqt files
do
     echo ${SYSTEM}
     #This creates a working directory for each system
     mkdir ${new_dir}/${SYSTEM}
done
cp Submit_Scripts/* ${new_dir}
