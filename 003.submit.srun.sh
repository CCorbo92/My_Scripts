#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=AutoDock
#SBATCH --output=AutoDock_dock
#SBATCH -p rn-long
system_file="/gpfs/projects/rizzo/ccorbo/DOCK6_with_ambpdb/SB_2020_testset/clean.systems.all "
new_dir="2021_06_13_srun"
for system in `cat ${system_file}`; do
   srun --exclusive -N1 -n1  bash ../run003.AutoDock4.docking.sh ${system} /gpfs/scratch/ccorbo/2021_05_autodock/AutoDock4_Tutorial ${new_dir}&> ./output/name.${system}.out    &
done

wait

