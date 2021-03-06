list_of_fam="all_DN_systems.txt"
for ref_fam in `cat ${list_of_fam}`; do
mkdir ${ref_fam}
cd ${ref_fam}
for i in {1..24}; do
mkdir anc_${i}
cd anc_${i}
cat >dn_generic.in<<EOF
conformer_search_type                                        denovo
dn_fraglib_scaffold_file                                     /gpfs/projects/rizzo/ccorbo/dock6.1_16_21_double_end_postprocess_new_beta/parameters/fraglib_scaffold.mol2
dn_fraglib_linker_file                                       /gpfs/projects/rizzo/ccorbo/dock6.1_16_21_double_end_postprocess_new_beta/parameters/fraglib_linker.mol2
dn_fraglib_sidechain_file                                    /gpfs/projects/rizzo/ccorbo/dock6.1_16_21_double_end_postprocess_new_beta/parameters/fraglib_sidechain.mol2
dn_user_specified_anchor                                     yes
dn_fraglib_anchor_file                                       /gpfs/scratch/ccorbo/MW_DN_tests/anchors/anchor_${i}.mol2
dn_torenv_table                                              /gpfs/projects/rizzo/ccorbo/dock6.1_16_21_double_end_postprocess_new_beta/parameters/fraglib_torenv.dat
dn_use_roulette                                              yes
dn_name_identifier                                           denovo
dn_sampling_method                                           graph
dn_graph_max_picks                                           30
dn_graph_breadth                                             3
dn_graph_depth                                               2
dn_graph_temperature                                         100.0
dn_pruning_conformer_score_cutoff                            100.0
dn_pruning_conformer_score_scaling_factor                    1.0
dn_pruning_clustering_cutoff                                 100.0
dn_upper_constraint_mol_wt                                   550
dn_lower_constraint_mol_wt                                   300
dn_constraint_mol_wt                                         550.0
dn_mol_wt_cutoff_type                                        soft
dn_mol_wt_std_dev                                            35
dn_constraint_rot_bon                                        15
dn_constraint_formal_charge                                  2.0
dn_heur_unmatched_num                                        1
dn_heur_matched_rmsd                                         2.0
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
internal_energy_cutoff                                       100.0
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           /gpfs/projects/rizzo/SB2012_testset_Built_2020/${ref_fam}/${ref_fam}.rec.clust.close.sph
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
grid_score_grid_prefix                                       /gpfs/projects/rizzo/SB2012_testset_Built_2020/${ref_fam}/${ref_fam}.rec
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
simplex_trans_step                                           0.20
simplex_rot_step                                             0.1
simplex_tors_step                                            15.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  250
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/flex_drive.tbl
EOF

cd ../
done
cd ../
done
