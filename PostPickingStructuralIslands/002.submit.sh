#!/bin/sh 
#SBATCH --partition=rn-long
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name=hmscalc
#SBATCH --output=testset.out
./002.pairwise_hms.sh
