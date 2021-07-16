#!/bin/sh
#SBATCH --partition=rn-long
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=CD_DOCK
#SBATCH --output=CD_DOCK.out

# This script calls RMSDh_extract.py which gets the RMSD for cartesian minimization and all output poses for each crossdocking pair

module unload anaconda/2
module load anaconda/3

DOCK_DIR="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking"
WORK_DIR="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/21_07_06"
CROSSDOCK_DIR="${WORK_DIR}/zzz.crossdock"
LIST_DIR="${DOCK_DIR}/zzz.sample_lists"
SCRIPT_DIR="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock"
cd ${CROSSDOCK_DIR}

echo "Getting results for family: "

list_of_fam="${LIST_DIR}/family_more_than7.txt"
for ref_fam in `cat ${list_of_fam}`; do  ### Open for loop 1
echo ${ref_fam} 
cd ${ref_fam}
rm *outcome.txt

list_of_sys1="${LIST_DIR}/${ref_fam}.txt"
for ref_system in `cat ${list_of_sys1}`; do
cd ${ref_system}

list_of_sys2="${LIST_DIR}/${ref_fam}.txt"
for comp_system in `cat ${list_of_sys2}`; do
cd ${comp_system}

python ${SCRIPT_DIR}/RMSDh_extract.py ${DOCK_DIR}/zzz.crossdock/${ref_fam}/${comp_system}/${ref_system}/${ref_system}_${comp_system}.min_scored.mol2 21_07_07_srun/summary_of_results_1.0 > ${comp_system}_${ref_system}.outcome.txt

echo  "$(<${comp_system}_${ref_system}.outcome.txt )"  >> ../../${ref_system}.outcome.txt

cd .. #back to outer system
done
cd .. # Back to family
done
cd .. #Back to base directory
done
