#!/bin/sh
#SBATCH --partition=rn-long
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=ADCD_Set4
#SBATCH --output=ADCD_Set4

# This script runs the cartesian minimization and docking set up in the last step. Make sure to specify the set being used in line below.
WORK_DIR="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/21_07_06/zzz.crossdock"
cd ${WORK_DIR}

list_of_fam="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/set_4.txt"
#ref_fam=`cat /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/family_more_than7.txt | head -n1`
for ref_fam in `cat ${list_of_fam}`; do
cd ${WORK_DIR}/${ref_fam}
echo -n "Running Family: "
echo ${ref_fam}

list_of_sys="/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/CrossDocking/zzz.sample_lists/${ref_fam}.txt"
for comp_system in `cat ${list_of_sys}`; do
cd ${WORK_DIR}/${ref_fam}/${comp_system}
#srun --exclusive -N1 -n1 -W 0  /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/run003.AutoDock4.docking.sh ${ref_fam}  ${comp_system} 21_07_07_srun &
srun --exclusive -N1 -n1  /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/CrossDocking_AutoDock/run003.AutoDock4.docking.sh ${ref_fam}  ${comp_system} 21_07_07_srun &

done

wait

done

