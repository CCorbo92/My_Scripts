#!/bin/sh
#SBATCH --partition=rn-long-40core
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=simpbuild
#SBATCH --output=simpbuild.out

for i in {1..20}; do
mkdir anc_${i}
cd anc_${i}
cat > dn_generic.in  << EOF
conformer_search_type                                        denovo
dn_fraglib_scaffold_file                                     /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Layered_FP/fraglib/fraglib_scaffold.mol2
dn_fraglib_linker_file                                       /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Layered_FP/fraglib/fraglib_linker.mol2
dn_fraglib_sidechain_file                                    /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Layered_FP/fraglib/fraglib_sidechain.mol2
dn_user_specified_anchor                                     yes
dn_fraglib_anchor_file                                       /gpfs/projects/rizzo/ccorbo/Testing_Grounds/srun_dn/anchors/anchor_${i}.mol2
dn_torenv_table                                              /gpfs/projects/rizzo/ccorbo/DOCK_Builds/Pushing/alpha_branch/dock6/parameters/fraglib_torenv.dat
dn_use_roulette                                              no
dn_name_identifier                                           denovo
dn_sampling_method                                           rand
dn_num_random_picks                                          25
dn_pruning_conformer_score_cutoff                            100
dn_pruning_conformer_score_scaling_factor                    1
dn_pruning_clustering_cutoff                                 100
dn_mol_wt_cutoff_type                                        hard
dn_upper_constraint_mol_wt                                   550.0
dn_lower_constraint_mol_wt                                   0
dn_constraint_rot_bon                                        15
dn_constraint_formal_charge                                  2
dn_heur_unmatched_num                                        1
dn_heur_matched_rmsd                                         2
dn_unique_anchors                                            1
dn_max_grow_layers                                           9
dn_max_root_size                                             25
dn_max_layer_size                                            25
dn_max_current_aps                                           5
dn_max_scaffolds_per_layer                                   1
dn_write_checkpoints                                         yes
dn_write_prune_dump                                          no
dn_write_orients                                             no
dn_write_growth_trees                                        no
dn_output_prefix                                             output
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
grid_score_primary                                           no
gist_score_primary                                           no
multigrid_score_primary                                      no
dock3.5_score_primary                                        no
continuous_score_primary                                     no
footprint_similarity_score_primary                           no
pharmacophore_score_primary                                  no
hbond_score_primary                                          no
interal_energy_score_primary                                 yes
internal_energy_rep_exp                                      12
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           0.25
simplex_rot_step                                             0.1
simplex_tors_step                                            15.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  250
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Layered_FP/21_12_18_AnchorTestingSecMin/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/rizzo/ccorbo/DOCK_Builds/Pushing/alpha_branch/dock6/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/rizzo/ccorbo/DOCK_Builds/Pushing/alpha_branch/dock6/parameters/flex_drive.tbl
EOF

/gpfs/projects/rizzo/ccorbo/DOCK_Builds/Pushing/PushToBetaBranch/compiled_chris_beta_branch/dock6/bin/dock6 -i dn_generic.in
cd ..
done
