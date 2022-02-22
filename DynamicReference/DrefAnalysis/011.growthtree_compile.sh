list_of_fam="smalltest10.txt"

for i in `seq 1 9`; do
for ref_sys in `cat ${list_of_fam}`; do
cd ${ref_sys}
rm All_Anchors_Layer_${i}.mol2

cat */output.anchor_1.root_layer_${i}.mol2 >>  All_Anchors_Layer_${i}.mol2

cat > hungarian.in << EOF
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             All_Anchors_Layer_${i}.mol2
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
descriptor_use_grid_score                                    yes
descriptor_use_pharmacophore_score                           no
descriptor_use_tanimoto                                      yes
descriptor_use_hungarian                                     yes
descriptor_use_volume_overlap                                yes
descriptor_use_gist                                          no
descriptor_use_dock3.5                                       no
descriptor_grid_score_rep_rad_scale                          1
descriptor_grid_score_vdw_scale                              1
descriptor_grid_score_es_scale                               1
descriptor_grid_score_grid_prefix                            /gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset/${ref_sys}/${ref_sys}.rec
descriptor_fingerprint_ref_filename                          /gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset/${ref_sys}/${ref_sys}.lig.am1bcc.mol2
descriptor_hms_score_ref_filename                            /gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset/${ref_sys}/${ref_sys}.lig.am1bcc.mol2
descriptor_hms_score_matching_coeff                          -5
descriptor_hms_score_rmsd_coeff                              1
descriptor_volume_score_reference_mol2_filename              /gpfs/projects/rizzo/ccorbo/Building_Stuff/DOCK6_with_ambpdb/SB_2020_testset/${ref_sys}/${ref_sys}.lig.am1bcc.mol2
descriptor_volume_score_overlap_compute_method               analytical
descriptor_weight_grid_score                                 0
descriptor_weight_fingerprint_tanimoto                       0
descriptor_weight_hms_score                                  1
descriptor_weight_volume_overlap_score                       0
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/rizzo/ccorbo/Testing_Grounds/Layered_FP/21_12_18_AnchorTestingSecMin/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3/parameters/flex_drive.tbl
chem_defn_file                                               /gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3/parameters/chem.defn
ligand_outfile_prefix                                        ${ref_sys}_rescored_fullref_Layer_${i}
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 yes
max_ranked_ligands                                           1
EOF

/gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6.21_11_05_Dynamic_Reference_Developmental_V1.5/bin/dock6 -i hungarian.in -o hungarian.out

echo "Finished ranking"

fragstring=`grep -A1 "USER_CHARGES" ${ref_sys}_rescored_fullref_Layer_${i}_ranked.mol2 | tail -n1`
echo ${fragstring}

file=`grep ${fragstring} */output.growth_tree*mol2 | awk -F ':' '{print $1}' | head -n1`
echo ${file}

cp ${file} ${ref_sys}_Best_Mol_Layer_${i}_GrowthTree.mol2

cd ..
done
done
