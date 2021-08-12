#!/bin/tcsh -fe

#
# This script should be executed after all standard dock to grid (run.006a) jobs have completed.
# It will first identify molecules that did dock - either successfully or unsuccessfully (Could
# not complete growth). Then it will take the diff between those molecules and the list of all
# original molecules. Undocked molecules (typically from jobs that hit the wall clock limit or
# that died spontaneously) will be concatenated into a new, "leftover" chunk. The user should
# specify the sizelim - or the minimum number of "leftover" molecules they would consider taking
# the time to dock. The numprocs variable should be set to the number of processors that were used
# on the dock to grid jobs.
#


### Set some variables manually
set sizelim  = "5000"
set numprocs = "1024"

### Set some paths
set dockdir   = "${DOCKHOME}/bin"
set amberdir  = "${AMBERHOME}/bin"
set moedir    = "${MOEHOME}/bin"
set rootdir   = "${VS_ROOTDIR}"
set mpidir    = "${VS_MPIDIR}"

set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"


### Choose parameters for cluster
### LIRED    24 ppn
### SeaWulf  28 ppn
### Rizzo    24 ppn
set wcl   = 48:00:00
set days = 0
set nodes = 4
set ppn   = 28
set queue = "long"
@ np      = (${nodes} * ${ppn})


### Make the appropriate directory
if (! -e ${rootdir}/${system}/006.dock-to-grid/${vendor}/) then
	echo "It is too soon to run this script! Finish the dock-to-grid first."
	exit
endif

rm -rf ${rootdir}/${system}/006.dock-to-grid/${vendor}/leftover/
mkdir -p ${rootdir}/${system}/006.dock-to-grid/${vendor}/leftover/
cd ${rootdir}/${system}/006.dock-to-grid/${VS_VENDOR}/leftover/


### Count the number of chunks
set num_chunks = `ls -l ${rootdir}/${system}/005.zinclibs/${vendor}/ | grep chunk | grep mol2 | wc -l`
set chunk = "0"
echo "num_chunks = ${num_chunks}"
echo "Counting and identifying undocked molecules..."


### Iterate over each chunk
while (${chunk} < ${num_chunks})


	### Identify all zinc ids in the original chunk 
	grep "Name" ${rootdir}/${system}/005.zinclibs/${vendor}/chunk${chunk}_scored.mol2 | awk '{print $3}' | sort | uniq > all_zincids.txt


	### Identify all zinc ids that were successfully docked
	grep "ZINC" ../chunk${chunk}/${system}.${vendor}.${chunk}.dock_to_grid.out | awk '{print $2}' > docked_zincids.txt


	### Identify the zinc ids that were docked, but growth could not be completed
	set counter = "1"
	while (${counter} < ${numprocs})
		grep -B 5 "Could" ../chunk${chunk}/${system}.${vendor}.${chunk}.dock_to_grid.out.${counter} | grep "ZINC" | awk '{print $2}' >> cncg_zincids.txt
		@ counter++
	end


	### Concatenate those last two chunks together and sort
	cat cncg_zincids.txt >> docked_zincids.txt
	cat docked_zincids.txt | sort | uniq > docked_zincids_final.txt


	### Find out which zinc ids were never docked
	diff all_zincids.txt docked_zincids_final.txt | grep "ZINC" | awk '{print $2}' > diff_zincids_chunk${chunk}.txt
	cat diff_zincids_chunk${chunk}.txt >> leftover.txt


	### If there are differences, add this chunk to the list of chunks
	if ( ` grep -c "ZINC" diff_zincids_chunk${chunk}.txt ` > 0) then
                echo "chunk${chunk}" >> chunk_list.txt
        endif

	rm all_zincids.txt docked_zincids.txt cncg_zincids.txt docked_zincids_final.txt
	@ chunk++
end


### Count the number of leftover molecules
set num_leftover = ` wc -l diff_zincids_chunk*.txt | grep "total" | awk '{print $1}' `
set chunk = "0"


### Check if the number is above the sizelim threshold
if (${num_leftover} < ${sizelim}) then

	echo "There are ${num_leftover} un-docked molecules, which is below the threshold of ${sizelim}"
else

	echo "There are ${num_leftover} un-docked molecules that will be resubmitted"


	### Break each chunk with undocked molecules into individual mol2s
	mkdir temp/
	cd temp/
	foreach multi_mol2 ( ` cat ../chunk_list.txt ` )
		python ${scriptdir}/break_into_mol.py ${rootdir}/${system}/005.zinclibs/${vendor}/${multi_mol2}_scored.mol2
	end

	foreach mol2 ( ` cat ../leftover.txt ` )
		cat ${mol2}.mol2 >> ../leftover.mol2
	end
	cd ../
	rm -r temp/
	
	echo "Docking chunk named leftover.mol2"


### Write the dock.in file
##################################################
cat <<EOF >${system}.${vendor}.leftover.dock_to_grid.in
conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               100.0
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            no
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             leftover.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ${rootdir}/${system}/003.spheres/${system}.rec.clust.close.sph
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ${rootdir}/${system}/004.grid/${system}.rec
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
ligand_outfile_prefix                                        ${vendor}.leftover.output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
##################################################


### Write the Cluster submit file for Maui
##################################################
if (`hostname -f` == "login1.cm.cluster"|| `hostname -f` == "login2.cm.cluster" ) then

cat <<EOF >${system}.${vendor}.leftover.qsub.sh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -q ${queue}
#PBS -N ${system}.${vendor}.leftover
#PBS -V

cd ${rootdir}/${system}/006.dock-to-grid/${vendor}/leftover/

${mpidir}/mpirun -np ${np} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.leftover.dock_to_grid.in \
-o ${system}.${vendor}.leftover.dock_to_grid.out

EOF

##################################################

### Write the Cluster submit file for slurm
##################################################
###(`hostname -f` == "rizzo.cm.cluster" )
else
cat <<EOF >${system}.${vendor}.leftover.qsub.sh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH -job-name=${system}.${vendor}.leftover
#SBATCH --output=${system}.${vendor}.leftover.out
#SBATCH -p ${queue}

echo "Job started on"
date
cd ${rootdir}/${system}/006.dock-to-grid/${vendor}/leftover/

${mpidir}/mpirun -np ${np} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.leftover.dock_to_grid.in \
-o ${system}.${vendor}.leftover.dock_to_grid.out

echo "Job finished on"
date

EOF

endif
##################################################



        ### Submit the job
        echo "Submitting ${system}.${vendor}.leftover.dock_to_grid " >> ../zzz.submit.log
        qsub ${system}.${vendor}.leftover.qsub.sh > & ${system}.${vendor}.leftover.qsub.log
        date >> ../zzz.submit.log

        cd ../
endif


