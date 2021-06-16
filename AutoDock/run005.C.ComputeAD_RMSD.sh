export system_file=$1
export docking_dir=$2

for system in `  cat ${system_file} `
do
echo ${system}

cd ${system}/${docking_dir}

rm RMSD_AD_grep.txt

grep " RANKING" ${system}.docking.dlg | awk '{print $6}' >> RMSD_AD_grep.txt


cd ../../

done
