#!/bin/tcsh -fe

#
#


### Set some variables manually 
set max_size = "500"
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
if (! -e ${rootdir}/${system}/011.final-results/) then
	mkdir -p ${rootdir}/${system}/011.final-results/
endif

if (! -e ${rootdir}/${system}/011.final-results/${vendor}) then
	mkdir -p ${rootdir}/${system}/011.final-results/${vendor}
endif

if (! -e ${rootdir}/${system}/011.final-results/system-files/) then
        echo "Creating the system-files directory and copying the corresponding files\n"
	mkdir -p ${rootdir}/${system}/011.final-results/system-files/
	cd ${rootdir}/${system}/011.final-results/system-files/
	cp ${rootdir}/${system}/001.lig-prep/${system}.lig.am1bcc.mol2 ./
	cp ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.python.min.mol2 ./
	cp ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2 ./
	cp ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.pdb ./
	cp ${rootdir}/${system}/002.rec-prep/pro.noH.pdb ./${system}.rec.noH.pdb
	cp ${rootdir}/${system}/003.spheres/${system}.rec.clust.close.sph ./
	cp ${rootdir}/${system}/004.grid/box.pdb ./${system}.box.pdb
endif

foreach dockscore (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: )

  echo ${dockscore}

rm -rf ${rootdir}/${system}/011.final-results/${vendor}/${dockscore}_rank
mkdir -p ${rootdir}/${system}/011.final-results/${vendor}/${dockscore}_rank
cd ${rootdir}/${system}/011.final-results/${vendor}/${dockscore}_rank
cat <<EOF> ${dockscore}.submit.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=output_submit1
#SBATCH --output=output_submit1
#SBATCH -p ${queue}





### Write out mol2 files for each different sorting method and for each group - families and 
### clusterheads. Each resulting mol2 file will contain $max_size molecules.
   echo " Writing the families and clusterhead mol2 files for each scoring metric\n"

foreach second_score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: )
	foreach group (clusterheads families)
	cp ${rootdir}/${system}/010.rdkit-and-cluster/${vendor}/${dockscore}_rank/\${group}_\${second_score}_secondary.csv ./
        python ${scriptdir}/rewrite_csv.py \${group}_\${second_score}_secondary.csv >> \${group}_\${second_score}_secondary_nocarr.csv
        set i="1"
	while (\$i < 1001)
        ###GET RID OF THE TRAILING .0!!!!!!!!!!!!!!!!!!!!!!!!!!!
        set linestart=\`sed -n "\${i}p" \${group}_\${second_score}_secondary_nocarr.csv | awk -F',' '{print \$50 }'| rev | cut -c 3- | rev \`
        set linestop=\`sed -n "\${i}p" \${group}_\${second_score}_secondary_nocarr.csv | awk -F',' '{print \$51 }'| rev | cut -c 3- | rev \`
        set chunk=\`sed -n "\${i}p" \${group}_\${second_score}_secondary_nocarr.csv | awk -F',' '{print \$52 }' \`
        set linequit=\${linestop}
        @ linequit++
        echo \${linestart} " " \${linestop} " "  \${linequit} " "  \${chunk} "  EOL">> \${group}_\${second_score}_XtrctInfo.txt
        sed -n "\${linestart},\${linestop}p;\${linequit}q"  ${rootdir}/${system}/008.rescore_mol/${vendor}/descriptor_score_rank/chunks/\${chunk}_ranked.mol2 >> \${group}_${dockscore}_\${second_score}_1000.mol2
	@ i++
        end
end
end
EOF
end


