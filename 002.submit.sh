#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=AutoDock
#SBATCH --output=AutoDock
#SBATCH -p rn-long

bash ../run002.AutoDock4.grid.generation.sh /gpfs/projects/rizzo/ccorbo/DOCK6_with_ambpdb/SB_2020_testset/clean.systems.all ./Tutorial.3.Score.csv
