#!/bin/sh 
#SBATCH --partition=long-24core
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name=Spike
#SBATCH --output=testset.out

./Hungarian_Cluster.sh
