#!/bin/sh 
#SBATCH --partition=rn-long
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name=desc_calc
#SBATCH --output=desc_calc.out
./002.pairwise_desc.sh
