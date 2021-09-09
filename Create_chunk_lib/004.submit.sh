#!/bin/sh 
#SBATCH --partition=rn-long-40core
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=ZINC_GET
#SBATCH --output=ZINC_GET.out

./run.004.make_chunks.sh
