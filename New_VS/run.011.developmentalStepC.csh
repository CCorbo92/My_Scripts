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
mkdir ${rootdir}/${system}/011.final-results/${vendor}/Final_Files
foreach dockscore (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: )
  echo ${dockscore}
mkdir ${rootdir}/${system}/011.final-results/${vendor}/Final_Files/${dockscore}
cd ${rootdir}/${system}/011.final-results/${vendor}/${dockscore}_rank

foreach second_score (Descriptor_Score: Continuous_Score: Pharmacophore_Score: Hungarian_Matching_Similarity_Score: Property_Volume_Score: desc_FPS_vdw_fps: desc_FPS_es_fps: Footprint_Similarity_Score: )
	foreach group (clusterheads families)
rm -r temp
mkdir temp
cd temp
csplit -q -z --digits=4 ../RDKit_FP.output_${group}.${second_score}_scored.mol2 /Name/ '{*}'
#for each file in the directory
ls >> list_of_files_tmp.txt
grep -v "list_of_files_tmp.txt" list_of_files_tmp.txt | grep -v "xx0000" | sort -n -r >> list_of_files.txt
cp ../${group}_${dockscore}_${second_score}_1000.mol2 ./
grep -n "Internal_energy_repulsive:" ${group}_${dockscore}_${second_score}_1000.mol2 | awk -F':' '{print $1}' >> insert_lines_nums_tmp.txt
#Do in reverse order
sort -n -r insert_lines_nums_tmp.txt > insert_lines_nums.txt
tac ../${group}_${second_score}_secondary.csv >> ./clusterheads_${second_score}_secondary_rev.csv

set i="1"
set length=`wc -l list_of_files.txt |awk '{print $1}'`
echo $length

@ length++
        while ($i < $length)
        set filename=`sed -n "${i}p" list_of_files.txt`
        echo ${filename}
        head -n 19 ${filename} | grep -v "DOCK" | grep -v "Molecular" | grep -v "Formal" > tmpheader.txt

        #add the cluster information here to the tmpheader file
        set clustrnum=`sed -n "${i}p"  ./clusterheads_${second_score}_secondary_rev.csv| awk -F ',' '{print $2}' `
        set clustrsiz=`sed -n "${i}p"  ./clusterheads_${second_score}_secondary_rev.csv| awk -F ',' '{print $3}' `
        echo "##########                         Cluster_Num:                   ${clustrnum}" >> tmpheader.txt
        echo "##########                        Cluster_Size:                   ${clustrsiz}" >> tmpheader.txt

        set linenum=`sed -n "${i}p" insert_lines_nums.txt`
        echo ${linenum}
        set molname=` echo ${group}_${dockscore}_${second_score}_1000.mol2`
        echo ${group}_${dockscore}_${second_score}_1000.mol2
        #insert tmpheader.txt into ${group}_${dockscore}_${second_score}_1000.mol2 at ${line}
        sed -e "${linenum}r tmpheader.txt" ${molname} > tmp.mol2
        mv tmp.mol2 ${group}_${dockscore}_${second_score}_1000.mol2
        @ i++
        end
mv ${group}_${dockscore}_${second_score}_1000.mol2 ../../Final_Files/${dockscore}/${group}_${dockscore}_${second_score}_1000.mol2
cd ..
pwd
rm -r temp
#end the for each
##########
end
end
end


