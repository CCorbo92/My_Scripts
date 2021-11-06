#!/bin/sh
#SBATCH --partition=rn-long-40core
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=35
#SBATCH --job-name=dnsrun
#SBATCH --output=dnsrun.out


for i in {1..35}; do
   srun --exclusive -N1 -n1 002.run_dn.sh  ${i} &
done

wait


