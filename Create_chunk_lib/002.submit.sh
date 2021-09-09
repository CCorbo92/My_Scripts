#!/bin/sh 
#SBATCH --partition=rn-long
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=Lib_Filter
#SBATCH --output=Lib_filter.out

./run.002.filter.sh 
