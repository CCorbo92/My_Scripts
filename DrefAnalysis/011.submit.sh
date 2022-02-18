#!/bin/sh
#SBATCH --partition=rn-long
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=rescore
#SBATCH --output=rescore.out


./011.growthtree_compile.sh
