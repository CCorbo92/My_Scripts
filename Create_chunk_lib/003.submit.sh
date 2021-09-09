#!/bin/sh 
#SBATCH --partition=rn-long
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=ZINC_GET
#SBATCH --output=ZINC_GET.out

./run.003.make_sets.sh
