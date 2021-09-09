#!/bin/tcsh 
#SBATCH --partition=rn-long-40core
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=ZINC_GET
#SBATCH --output=ZINC_GET.out

tcsh run.001.getfiles_from_zinc15.csh 
