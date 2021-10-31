#!/bin/sh 
#SBATCH --partition=rn-long-40core
#SBATCH --time=130:00:00
#SBATCH --nodes=3
#SBATCH --ntasks=40
#SBATCH --job-name=LigE3KL6
#SBATCH --output=LigEnrch.out

# This script will dock actives and decoys from the downloaded set of actives and decoys from DUDE
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 10/2021
# Last Edit by: Christopher Corbo

#Set the three variables below
testset="/gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset"
processes="120"
system="3KL6"
system_file="${system}_set.txt"
for systemrec in `cat ${system_file}`; do
./FLX_actives.sh ${system} ${systemrec} ${testset} ${processes}
./FLX_decoys.sh ${system} ${systemrec} ${testset} ${processes}

echo ${system} " has finished processing" >> DUDE_Docked_Log.txt
done

