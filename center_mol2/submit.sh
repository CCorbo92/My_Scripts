#!/bin/sh
#SBATCH --partition=rn-long
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --job-name=qiktst
#SBATCH --output=qiktest.out


module load chimera/1.13.1
mol_count="$(grep 'MOLECULE' ${1} | wc -l)"

for i in $(seq 1 $mol_count); do
   srun --exclusive -N1 -n1 center_srun.sh ${1} ${i} &
done

wait

sleep 20s

for i in $(seq 1 $mol_count); do
   cat ${i}_centered.mol2 >> centered_${1}
   rm ${i}_centered.mol2
done

grep -v "Step created for job" ./qiktest.out | grep -v "step creation temporarily disabled" | grep "cat:" | awk '{print $2}' |tr -d '_centered.mol2:' >> leftover.txt

leftover_file="leftover.txt"

for i in `cat ${leftover_file}`; do
      srun --exclusive -N1 -n1 center_srun.sh ${1} ${i} &
done

wait

sleep 20s


for i in `cat ${leftover_file}`; do
   cat ${i}_centered.mol2 >> centered_${1}
   rm ${i}_centered.mol2
done
