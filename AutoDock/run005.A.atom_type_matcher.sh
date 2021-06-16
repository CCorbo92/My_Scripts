export system_file=$1
export docking_dir=$2

for system in `  cat ${system_file} `
do
echo ${system}

cd ${system}/${docking_dir}
rm ${system}.docking.noH.modified.mol2
python  ../../../atom_type_matcher.py ${system} >> ${system}.docking.noH.modified.mol2

cd ../../

done
