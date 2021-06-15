export system_file=$1
export docking_dir=$2

for system in `  cat ${system_file} `
do
echo ${system}

cd ${system}/${docking_dir}

cat >rescore.in<<EOF
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ${system}.docking.noH.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes
rmsd_reference_filename                                      ${system}.lig.gast.noH.mol2
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              no
atom_model                                                   united
vdw_defn_file                                                /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/AMS536/zzz.programs/dock6.9_release/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ${system}.lig.AD
write_orientations                                           no
num_scored_conformers                                        10
cluster_conformations                                        no
rank_ligands                                                 no
EOF

dock6 -i rescore.in -o rescore.out

grep "HA_RMSDh:" ${system}.lig.AD_scored.mol2 | awk '{print $3}' >> RMSD_grep.txt

cd ../../

done
