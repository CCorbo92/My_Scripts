list_of_fam="../smalltest10.txt"
cd Histogram_Data_Dyn
for ref_sys in `cat ${list_of_fam}`; do

#Get best score on hms and average of 100 best from all anchors
grep "desc_HMS_num_Ref_hvy_atms_Mtchd" ${ref_sys}_rescored_fullref_scored.mol2 | awk '{print $3}' | sort -n >> ${ref_sys}_Ref_Match.txt

grep "desc_HMS_num_Pose_hvy_atms_Mtchd" ${ref_sys}_rescored_fullref_scored.mol2 | awk '{print $3}' | sort -n >> ${ref_sys}_Pose_Match.txt

grep "desc_HMS_num_Pose_hvy_atms_Unmtchd" ${ref_sys}_rescored_fullref_scored.mol2 | awk '{print $3}' | sort -n >> ${ref_sys}_Pose_Unmatch.txt

grep  "desc_HMS_rmsd_mtchd_hvy_atms" ${ref_sys}_rescored_fullref_scored.mol2 | awk '{print $3}' | sort -n >> ${ref_sys}_RMSD_Match.txt

done
