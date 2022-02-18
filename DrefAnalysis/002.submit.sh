#!/bin/sh
#SBATCH --partition=rn-long-40core
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=dnsrun2
#SBATCH --output=dnsrun2.out


for i in {1..35}; do
   echo ${i}
   srun --exclusive -N1 -n1 002.run_dn.sh  ${i} &
done

wait


