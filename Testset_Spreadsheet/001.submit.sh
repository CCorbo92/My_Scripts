#!/bin/sh 
#SBATCH --partition=rn-long-40core
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Utils
#SBATCH --output=001Util.out

# This script will call 001.utils.calculate.sh which runs DOCK to calculate DOCK and RDKit descriptors

# Path to prepared testset files
testset="/gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset"
cd ${testset}
bash ./data_used_for_spreadsheet/001.utils.calculate.sh

