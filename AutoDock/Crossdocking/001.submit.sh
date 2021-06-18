#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=AutoDock
#SBATCH --output=AutoDock
#SBATCH -p rn-long

bash ../run001.AutoDock4.system.prep.sh ../clean.systems.all ./Tutorial.3.Score.csv
