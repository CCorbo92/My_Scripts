#!/bin/tcsh -fe

#
#


### Set some variables manually 
set max_size = "1000"
set max_num  = "100000"
set cutoff   = "0.2"
set max_res  = "50"


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

set wcl   = 7-00:00:00
set nodes = 1
set ppn   = 24
set queue = "rn-long"
### Make a directory for compiling all of the docked results for a given vendor. If the top
### directory already exists, don't remove other vendor results.

foreach dockscore (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: )

  echo ${dockscore}

cd ${rootdir}/${system}/011.final-results/${vendor}/${dockscore}_rank
cat <<EOF> ${dockscore}.RDkit.submit.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=${dockscore}.rdkit
#SBATCH --output=output_submit1
#SBATCH -p ${queue}




foreach second_score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: )
	foreach group (clusterheads families)
        # Write the RDkit dock input files here

cat <<EOF >\${group}.\${second_score}.RDKit_FP_dock.in
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             \${group}_${dockscore}_\${second_score}_1000.mol2
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
dbfilter_max_molwt                                           9999
dbfilter_min_molwt                                           0
dbfilter_max_formal_charge                                   10
dbfilter_min_formal_charge                                   -10
dbfilter_max_stereocenters                                   100
dbfilter_min_stereocenters                                   0
dbfilter_max_spiro_centers                                   100
dbfilter_min_spiro_centers                                   0
dbfilter_max_clogp                                           20
dbfilter_min_clogp                                           -20
dbfilter_max_logs                                            20
dbfilter_min_logs                                            -20
dbfilter_max_qed                                             1
dbfilter_min_qed                                             0
dbfilter_max_sa                                              10
dbfilter_min_sa                                              0
dbfilter_max_pns                                             100
filter_sa_fraglib_path                                       /gpfs/projects/rizzo/spak/builds/dock6_copy/parameters/sa_fraglib.dat
filter_PAINS_path                                            /gpfs/projects/rizzo/spak/builds/dock6_copy/parameters/pains_table.dat
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
grid_score_primary                                           no
gist_score_primary                                           no
multigrid_score_primary                                      no
dock3.5_score_primary                                        no
continuous_score_primary                                     no
footprint_similarity_score_primary                           yes
fps_score_use_footprint_reference_mol2                       yes
fps_score_footprint_reference_mol2_filename                  ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.python.min.mol2
fps_score_foot_compare_type                                  Euclidean
fps_score_normalize_foot                                     no
fps_score_foot_comp_all_residue                              yes
fps_score_receptor_filename                                  ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2
fps_score_vdw_att_exp                                        6
fps_score_vdw_rep_exp                                        12
fps_score_vdw_rep_rad_scale                                  1
fps_score_use_distance_dependent_dielectric                  yes
fps_score_dielectric                                         4
fps_score_vdw_fp_scale                                       1
fps_score_es_fp_scale                                        1
fps_score_hb_fp_scale                                        0
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/rizzo/spak/builds/dock6_copy/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/rizzo/spak/builds/dock6_copy/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/rizzo/spak/builds/dock6_copy/parameters/flex_drive.tbl
ligand_outfile_prefix                                        RDKit_FP.output_\${group}.\${second_score}
write_footprints                                             yes
write_hbonds                                                 no
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
\EOF
/gpfs/projects/rizzo/spak/builds/dock6_copy/bin/dock6 -i \${group}.\${second_score}.RDKit_FP_dock.in -o \${group}.\${second_score}.RDKit_FP_dock.out
end
end
EOF
sed -i 's/\\EOF/EOF/g' ${dockscore}.RDkit.submit.csh
end


