

cd ./mol2_decompress
ls *.mol2 > mol2_list.txt

mkdir ./filtered_mol2
cd ./filtered_mol2

### Here, begin writing a tcsh script that can be run on the queue
##################################################
count="0"
for i in `cat ../mol2_list.txt`; do
##############################################
cat > ${count}.dock_filter.in  << EOF  
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             ../$i
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          yes
dbfilter_max_heavy_atoms                                     9999
dbfilter_min_heavy_atoms                                     0
dbfilter_max_rot_bonds                                       15
dbfilter_min_rot_bonds                                       0
dbfilter_max_molwt                                           550.0
dbfilter_min_molwt                                           200.0
dbfilter_max_formal_charge                                   2
dbfilter_min_formal_charge                                   -2
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/rizzo/zzz.programs/dock6.9_release/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/rizzo/zzz.programs/dock6.9_release/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/rizzo/zzz.programs/dock6.9_release/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ${count}.output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
############################################################

/gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6_with_rdkit/bin/dock6 -i ${count}.dock_filter.in -o ${count}.dock_filter.out

count=$((count+1))

done
#################################################







