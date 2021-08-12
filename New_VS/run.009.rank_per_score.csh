#!/bin/tcsh -fe
### Steven please for the love of god make the footprint text file
### Add in the calculation of RDKit descriptors at this step
### Set some variables manually
set attractive   = "6" 
set repulsive    = "12"

set max_num = "${MAX_NUM_MOL}"

### Set some paths
set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set mpidir    = "${VS_MPIDIR}/bin"

set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set descriptdir = "${rootdir}/zzz.descriptor"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"


set wcl   = 7-00:00:00
set nodes = 1 
set ppn   = 40
set queue = "rn-long-40core"
@ numprocs = (${nodes} * ${ppn})


if ( ! -e ${rootdir}/${system}/008.rescore_mol/${vendor}/descriptor_score_rank/chunks ) then
        echo "Cartesian minimized and rescored mol2 chunks do not seem to exist. Exiting.";
        exit
endif

### Make the appropriate directory. If it already exists, remove previous dock results from only### the same vendor.
if (! -e ${rootdir}/${system}/009.rank_per_score) then
        mkdir -p ${rootdir}/${system}/009.rank_per_score/
endif

rm -rf ${rootdir}/${system}/009.rank_per_score/${vendor}
mkdir -p ${rootdir}/${system}/009.rank_per_score/${vendor}
cd ${rootdir}/${system}/009.rank_per_score/${vendor}


################################
#########1st PHASE##############
################################
cat <<EOF> ${system}.${vendor}.submit1.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=output_submit1
#SBATCH --output=output_submit1
#SBATCH -p ${queue}


#module load shared
module unload anaconda/2
module load anaconda/3

squeue -u stpak -t RUNNING
module list


python ${scriptdir}/phase1.py ${rootdir}/${system}/008.rescore_mol/${vendor}/descriptor_score_rank/chunks
EOF


################################
#########2nd PHASE##############
################################
cat <<EOF> ${system}.${vendor}.submit2.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=output_submit2
#SBATCH --output=output_submit2
#SBATCH -p ${queue}

#module load shared
module unload anaconda/2
module load anaconda/3
squeue -u stpak -t RUNNING
module list

python ${scriptdir}/phase2.py ${max_num}

foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)

    rm \${score}_sorted.csv 
end

EOF

################################
#########3rd PHASE##############
################################

foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)

cat <<EOF >${score}.csh
#!/bin/tcsh 

#module load shared
module unload anaconda/2
module load anaconda/3

squeue -u stpak -t RUNNING

module list
python ${scriptdir}/phase3.py ${rootdir}/${system}/008.rescore_mol/${vendor}/descriptor_score_rank/chunks ${score} ${max_num}

EOF
end

foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)

cat <<EOF >${score}.RDKIT.csh
#!/bin/tcsh 

#module load shared
module unload anaconda/3
module load anaconda/2

squeue -u stpak -t RUNNING

module list

/gpfs/projects/rizzo/spak/builds/dock6_copy/bin/dock6 -i ${score}.RDKit_FP_dock.in -o ${score}.RDKit_FP_dock.out

EOF
end


cat <<EOF>${system}.${vendor}.submit3.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=output_submit3
#SBATCH --output=output_submit3
#SBATCH -p ${queue}

#module load shared
module unload anaconda/2
module load anaconda/3
squeue -u stpak -t RUNNING
module list

foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)

    srun --exclusive -N1 -n1 -W 0 \${score}.csh & 

end

wait
EOF



################################
#########4th PHASE##############
################################
foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)

cat <<EOF >${score}.dock.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ${score}_top_${max_num}.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          yes
dbfilter_max_heavy_atoms                                     9999
dbfilter_min_heavy_atoms                                     0
dbfilter_max_rot_bonds                                       9999
dbfilter_min_rot_bonds                                       0
dbfilter_max_hb_donors                                       999
dbfilter_min_hb_donors                                       0
dbfilter_max_hb_acceptors                                    999
dbfilter_min_hb_acceptors                                    0
dbfilter_max_molwt                                           9999.0
dbfilter_min_molwt                                           0.0
dbfilter_max_formal_charge                                   10.0
dbfilter_min_formal_charge                                   -10.0
dbfilter_max_stereocenters                                   6
dbfilter_min_stereocenters                                   0
dbfilter_max_spiro_centers                                   6
dbfilter_min_spiro_centers                                   0
dbfilter_max_clogp                                           20.0
dbfilter_min_clogp                                           -20.0
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                /gpfs/projects/rizzo/spak/builds/dock6_rdkit/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /gpfs/projects/rizzo/spak/builds/dock6_rdkit/parameters/flex.defn
flex_drive_file                                              /gpfs/projects/rizzo/spak/builds/dock6_rdkit/parameters/flex_drive.tbl
ligand_outfile_prefix                                        MACCS_${score}
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

cat <<EOF >${score}.RDKit_FP_dock.in
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             ${score}_top_${max_num}.mol2
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
ligand_outfile_prefix                                        RDKit_FP.output_${score}
write_footprints                                             yes
write_hbonds                                                 no
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF


end

cat <<EOF> ${system}.${vendor}.submit4.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=output_submit4
#SBATCH --output=output_submit4
#SBATCH -p ${queue}


#module load shared
squeue -u stpak -t RUNNING
module list

module unload anaconda/3
module load anaconda/2

foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)

     /gpfs/projects/rizzo/spak/builds/dock6_MACCS/bin/dock6 -i \${score}.dock.in -o \${score}.dock.out

     rm MACCS_\${score}*
     
     grep "MACCS:" \${score}.dock.out > MACCSraw_\${score}.txt
     rm \${score}.dock.in
     
     awk '{print \$2}' MACCSraw_\${score}.txt > MACCSfinal_\${score}.txt
     rm MACCSraw_\${score}.txt
     
     grep "Molecule:" \${score}.dock.out > ZINCraw_\${score}.txt
     awk '{print \$2}' ZINCraw_\${score}.txt > ZINCfinal_\${score}.txt
     
     rm ZINCraw_\${score}.txt
     
     paste ZINCfinal_\${score}.txt MACCSfinal_\${score}.txt -d "," > ZINC_MACCS_\${score}.txt
     
     rm MACCSfinal_\${score}.txt
     rm ZINCfinal_\${score}.txt
     rm \${score}.dock.out
end

EOF

cat <<EOF> ${system}.${vendor}.submit5.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=output_submit5
#SBATCH --output=output_submit5
#SBATCH -p ${queue}


#module load shared
squeue -u stpak -t RUNNING
module list

module unload anaconda/3
module load anaconda/2


foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)

    srun --exclusive -N1 -n1 -W 0 \${score}.RDKIT.csh &

end
wait


rm *RDKit_FP_dock.in
rm *RDKit_FP_dock.out 

foreach score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: Total_Score:)
    rm RDKit_FP.output_\${score}_scored.mol2
end
EOF
################################
#######INITIATION PHASE#########
################################

cat <<EOF> ${system}.${vendor}.submit.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=output_submit
#SBATCH --output=output_submit
#SBATCH -p ${queue}


#module load shared
module unload anaconda/2
module load anaconda/3

squeue -u stpak -t RUNNING
module list

chmod +x *.csh
srun --exclusive -N1 -n1 -W 0 ${system}.${vendor}.submit1.csh
srun --exclusive -N1 -n1 -W 0 ${system}.${vendor}.submit2.csh 
srun --exclusive -N1 -n1 -W 0 ${system}.${vendor}.submit3.csh  
srun --exclusive -N1 -n1 -W 0 ${system}.${vendor}.submit4.csh
srun --exclusive -N1 -n1 -W 0 ${system}.${vendor}.submit5.csh

date
EOF

sbatch ${system}.${vendor}.submit.csh 
