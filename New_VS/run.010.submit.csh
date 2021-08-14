#!/bin/tcsh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=run.010
#SBATCH --output=run.010.out
#SBATCH -p rn-long-40core

tcsh run.010.snippet_test.csh
#tcsh run.010.unique_and_cluster.csh
