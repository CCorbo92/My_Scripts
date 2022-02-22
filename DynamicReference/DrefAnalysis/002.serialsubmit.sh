#!/bin/sh
#SBATCH --partition=rn-long
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=dnsrun
#SBATCH --output=dnsrun.out


for i in {1..35}; do
   echo ${i}
   ./002.run_dn.sh  ${i} 
done

wait


