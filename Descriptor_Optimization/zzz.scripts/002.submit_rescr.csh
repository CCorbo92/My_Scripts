#!/bin/tcsh
#SBATCH --partition=rn-long-40core
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=rscor
#SBATCH --output=rscor.out

set system_in = "1B9V"
set scripts_dir = "/gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/zzz.scripts"

cd /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/${system_in}

source ${scripts_dir}/run.000.set_env_vars.csh ${system_in}

tcsh ${scripts_dir}/Rescore.csh Actives
tcsh ${scripts_dir}/Rescore.csh Decoys

rm *out*
