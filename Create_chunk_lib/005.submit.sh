#!/bin/sh 
#SBATCH --partition=rn-long
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=ZINC_GET
#SBATCH --output=Chunks.out

./run.005.really_make_them.sh
