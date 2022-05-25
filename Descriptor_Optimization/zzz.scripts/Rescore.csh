module load anaconda/2

### Set some variables manually
set attractive   = "6"
set repulsive    = "9"
set grid_spacing = "0.3"
set box_margin   = "8"
set sphcut       = "8"
set maxkeep      = "75"


### Set some paths
set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"
set mpidir    = "${VS_MPIDIR}/bin"

set class = "${1}"

cat <<EOF >${system}.${class}_rescore.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/${system}/${class}_min_scored.mol2 
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     no
continuous_score_secondary                                   no
footprint_similarity_score_primary                           no
footprint_similarity_score_secondary                         no
pharmacophore_score_primary                                  no
pharmacophore_score_secondary                                no
descriptor_score_primary                                     yes
descriptor_score_secondary                                   no
descriptor_use_grid_score                                    no
descriptor_use_multigrid_score                               no
descriptor_use_continuous_score                              yes
descriptor_use_footprint_similarity                          yes
descriptor_use_pharmacophore_score                           yes
descriptor_use_tanimoto                                      yes
descriptor_use_hungarian                                     yes
descriptor_use_volume_overlap                                yes
descriptor_cont_score_rec_filename                           /gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset/${system}/${system}.rec.clean.mol2
descriptor_cont_score_att_exp                                6
descriptor_cont_score_rep_exp                                12
descriptor_cont_score_rep_rad_scale                          1
descriptor_cont_score_use_dist_dep_dielectric                yes
descriptor_cont_score_dielectric                             4.0
descriptor_cont_score_vdw_scale                              1
descriptor_cont_score_es_scale                               1
descriptor_fps_score_use_footprint_reference_mol2            yes
descriptor_fps_score_footprint_reference_mol2_filename       /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/${system}/Ref_min_scored.mol2
descriptor_fps_score_foot_compare_type                       Euclidean
descriptor_fps_score_normalize_foot                          no
descriptor_fps_score_foot_comp_all_residue                   yes
descriptor_fps_score_receptor_filename                       /gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset/${system}/${system}.rec.clean.mol2
descriptor_fps_score_vdw_att_exp                             6
descriptor_fps_score_vdw_rep_exp                             12
descriptor_fps_score_vdw_rep_rad_scale                       1
descriptor_fps_score_use_distance_dependent_dielectric       yes
descriptor_fps_score_dielectric                              4.0
descriptor_fps_score_vdw_fp_scale                            1
descriptor_fps_score_es_fp_scale                             1
descriptor_fps_score_hb_fp_scale                             0
descriptor_fms_score_use_ref_mol2                            yes
descriptor_fms_score_ref_mol2_filename                       /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/${system}/Ref_min_scored.mol2
descriptor_fms_score_write_reference_pharmacophore_mol2      no
descriptor_fms_score_write_reference_pharmacophore_txt       no
descriptor_fms_score_write_candidate_pharmacophore           no
descriptor_fms_score_write_matched_pharmacophore             no
descriptor_fms_score_compare_type                            overlap
descriptor_fms_score_full_match                              yes
descriptor_fms_score_match_rate_weight                       5.0
descriptor_fms_score_match_dist_cutoff                       1.0
descriptor_fms_score_match_proj_cutoff                       0.7071
descriptor_fms_score_max_score                               20
descriptor_fingerprint_ref_filename                          /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/${system}/Ref_min_scored.mol2
descriptor_hms_score_ref_filename                            /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/${system}/Ref_min_scored.mol2
descriptor_hms_score_matching_coeff                          -5
descriptor_hms_score_rmsd_coeff                              1
descriptor_volume_score_reference_mol2_filename              /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Benchmarking_and_Validation/LigandEnrichment/OptimizedDescriptor/${system}/Ref_min_scored.mol2
descriptor_volume_score_overlap_compute_method               analytical
descriptor_weight_cont_score                                 1
descriptor_weight_fps_score                                  2
descriptor_weight_pharmacophore_score                        4
descriptor_weight_fingerprint_tanimoto                       -80
descriptor_weight_hms_score                                  3
descriptor_weight_volume_overlap_score                       -80
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/rizzo/for_Pak_Chris/DOCK6_Screening_Protocols/zzz.parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/rizzo/for_Pak_Chris/DOCK6_Screening_Protocols/zzz.parameters/flex.defn
flex_drive_file                                              /gpfs/projects/rizzo/for_Pak_Chris/DOCK6_Screening_Protocols/zzz.parameters/flex_drive.tbl
chem_defn_file                                               /gpfs/projects/rizzo/for_Pak_Chris/DOCK6_Screening_Protocols/zzz.parameters/chem.defn
pharmacophore_defn_file                                      /gpfs/projects/rizzo/for_Pak_Chris/DOCK6_Screening_Protocols/zzz.parameters/ph4.defn
ligand_outfile_prefix                                        ${system}_${class}_rescore
write_footprints                                             no
write_hbonds                                                 no
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 yes
max_ranked_ligands                                           250000
EOF
##################################################

### Execute dock on the headnode
${mpidir}/mpirun -np 40 ${dockdir}/dock6.mpi -i ${system}.${class}_rescore.in -o ${system}.${class}_rescore.out
