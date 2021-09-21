#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=AutoDock
#SBATCH --output=AutoDock_dock
#SBATCH -p rn-long-40core
system_file=" List_all_seed_dir.txt"
for system in `cat ${system_file}`; do
  ./run005.C.ComputeAD_RMSD.sh ./clean.systems.all ${system}
  ./run006.calculateStats_AD.sh ./clean.systems.all ${system}
  grep "Success" Outcome_AD.txt | awk '{print $2}' >> Success_${system}.txt
  grep "Sco" Outcome_AD.txt | awk '{print $3}' >> Score_Fail_${system}.txt
  grep "Sam" Outcome_AD.txt | awk '{print $3}' >> Sample_Fail_${system}.txt
  mv Outcome_AD.txt Outcome_AD_${system}.txt

  ./run008.all_stats.sh ./clean.systems.all 2021_09_09_S4_4
  awk '{print $2}' timing_all.txt >> timing_all_only_${system}.txt
  mv timing_all.txt timing_all_${system}.txt
done
