#!/bin/sh 
#SBATCH --partition=rn-long-40core
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=Uniprot
#SBATCH --output=Uniprot.out
./001.access_page.sh
