#!/bin/sh
#SBATCH --partition=rn-long
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=ADCD
#SBATCH --output=ADCD

# This script runs the cartesian minimization and docking set up in the last step. Make sure to specify the set being used in line below.
WORK_DIR="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/21_07_06/zzz.crossdock"
cd ${WORK_DIR}

list_of_fam="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/family_more_than7.txt"
#ref_fam=`cat /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/family_more_than7.txt | head -n1`
for ref_fam in `cat ${list_of_fam}`; do
mkdir ${WORK_DIR}/${ref_fam}
cd ${WORK_DIR}/${ref_fam}
echo -n "Running Family: "
echo ${ref_fam}

list_of_sys="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/${ref_fam}.txt"
for comp_system in `cat ${list_of_sys}`; do
mkdir ${WORK_DIR}/${ref_fam}/${comp_system} 
cd ${WORK_DIR}/${ref_fam}/${comp_system}
srun --exclusive -N1 -n1 -W 0  /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/run002.AutoDock4.grid.generation.sh ${ref_fam}  ${comp_system} &

done

wait


done

