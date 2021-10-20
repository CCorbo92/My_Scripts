#!/usr/bin/sh


list_of_fam="viable.txt"
list_of_compfam="viable.txt"

for ref in `cat ${list_of_fam}`; do  ### Open for loop 1

mkdir ${ref}

for comp in `cat ${list_of_compfam}`; do
cd ${ref}

cat >rescore.in<<EOF
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ../Viable/${ref}
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
descriptor_use_continuous_score                              no
descriptor_use_footprint_similarity                          no
descriptor_use_pharmacophore_score                           no
descriptor_use_tanimoto                                      no
descriptor_use_hungarian                                     yes
descriptor_hms_score_ref_filename                            ../Viable/${comp}
descriptor_hms_score_matching_coeff                          -5
descriptor_hms_score_rmsd_coeff                              1
descriptor_weight_hms_score                                  1
descriptor_use_volume_overlap                                no
descriptor_fingerprint_ref_filename                          ../Viable/${comp}
descriptor_weight_fingerprint_tanimoto                       -5
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ${ref}_${comp}.output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

dock6 -i rescore.in 

hms=`grep  "Hungarian" ${ref}_${comp}.output_scored.mol2 |awk '{print $3}'`
val=`echo $hms'< -2.1' | bc -l `
if [ "$val" -eq "1" ]; then

name1=`grep "Name" ../Viable/${comp} | awk '{print $3}'`
name2=`grep "Name" ../Viable/${ref} | awk '{print $3}'`
echo ${name1} " " ${name2} >> ../PreclusterList.txt

fi

cd ..
done
done
