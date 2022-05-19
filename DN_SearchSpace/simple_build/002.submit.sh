#!/bin/sh
#SBATCH --partition=rn-long-40core
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=simpbuild
#SBATCH --output=simpbuild.out
for i in {1..20}; do
   echo ${i}
   srun --exclusive -N1 -n1 ./001.run_dn.sh  ${i} &
done
