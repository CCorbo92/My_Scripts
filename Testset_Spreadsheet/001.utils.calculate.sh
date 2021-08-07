# This script writes the dock input file to calculate descriptors with DOCK and RDKit 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 08/2021
# Last Edit by: Christopher Corbo

system_file="clean.systems.all"

for system in `cat ${system_file}`; do
cd ${system}

echo $system

cat >util.in<<EOF
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ${system}.lig.am1bcc.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          yes
dbfilter_max_heavy_atoms                                     999
dbfilter_min_heavy_atoms                                     0
dbfilter_max_rot_bonds                                       999
dbfilter_min_rot_bonds                                       0
dbfilter_max_hb_donors                                       999
dbfilter_min_hb_donors                                       0
dbfilter_max_hb_acceptors                                    999
dbfilter_min_hb_acceptors                                    0
dbfilter_max_molwt                                           9999.0
dbfilter_min_molwt                                           0.0
dbfilter_max_formal_charge                                   20.0
dbfilter_min_formal_charge                                   -20.0
dbfilter_max_stereocenters                                   999
dbfilter_min_stereocenters                                   0
dbfilter_max_spiro_centers                                   999
dbfilter_min_spiro_centers                                   0
dbfilter_max_clogp                                           40.0
dbfilter_min_clogp                                           -40.0
filter_sa_fraglib_path                                       /gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6_with_rdkit/parameters/sa_fraglib.dat
filter_PAINS_path                                            /gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6_with_rdkit/parameters/pains_table.dat
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6_with_rdkit/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6_with_rdkit/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6_with_rdkit/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ${system}.lig.descriptors
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

/gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6_with_rdkit/bin/dock6 -i util.in -o util.out

cd ..

done
