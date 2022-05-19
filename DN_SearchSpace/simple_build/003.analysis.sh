#!/bin/sh
#SBATCH --partition=rn-long-40core
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=rdkitscore
#SBATCH --output=rdkitscore.out
for i in {11..15}; do
cd anc_${i}

grep "\-\-" output.denovo_build.mol2 | grep -o -n '\.'  | cut -d : -f 1 | uniq -c| awk '{print $1}' >fragcount.txt
grep "Weight" output.denovo_build.mol2 | awk '{print $3}' > MW.txt

csplit -z output.denovo_build.mol2 /TEMP/ '{*}'
for file in xx*;
do

fragcount=`grep "Frag_String" $file | grep -o -n '\.'  | cut -d : -f 1 | uniq -c| awk '{print $1}'`
echo $fragcount
if [ $fragcount -gt 6 ]; then
   cat $file >> 7_up_frag.mol2
fi
done

rm xx*


cat >util.in<<EOF
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             7_up_frag.mol2
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

csplit -z util.out /Molecule:/ '{*}'
rm xx00

for out in xx*;
do
SA=`grep -w "SA" $out | awk '{print $2}' `
QED=`grep " QED" $out | awk '{print $2}' `
echo $SA
echo $QED

if (( $(echo "$QED > 0.49" |bc -l) )); then
if (( $(echo "$SA < 3.85" |bc -l) )); then
 cat $out >> passes_QED_SA.mol2
fi
fi

done

grep "TEMP"  7_up_frag.mol2 | wc -l >> ../Total_7_up_Built.txt
grep "Molecule:"  passes_QED_SA.mol2 | wc -l >> ../Total_Passing_QED_SA.txt

cd ..
done
