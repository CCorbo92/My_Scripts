#!/bin/sh
#SBATCH --partition=rn-long-40core
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=centering
#SBATCH --output=centering.out


module load chimera/1.13.1
system_file="clean.systems.all"

#Pass the command line argument of the mol2 you wish to center
mol_count="$(grep 'MOLECULE' All_testset.mol2 | wc -l)"

for i in $(seq 1 $mol_count); do
   srun --exclusive -N1 -n1 007.center_mol2.sh All_testset.mol2 ${i} &
done

wait

sleep 20s

for i in $(seq 1 $mol_count); do
   cat ${i}_centered.mol2 >> centered_All_testset.mol2
   rm ${i}_centered.mol2
done

grep -v "Step created for job" ./centering.out | grep -v "step creation temporarily disabled" | grep "cat:" | awk '{print $2}' |tr -d '_centered.mol2:' >> leftover.txt

leftover_file="leftover.txt"

for i in `cat ${leftover_file}`; do
      srun --exclusive -N1 -n1 007.center_mol2.sh All_testset.mol2 ${i} &
done

wait

sleep 20s


for i in `cat ${leftover_file}`; do
   cat ${i}_centered.mol2 >> centered_All_testset.mol2
   rm ${i}_centered.mol2
done

rm *out
