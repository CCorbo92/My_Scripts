source activate openbabel

export SYS=$1
export docking_dir=$2

for system in `  cat ${SYS} `
do
echo ${system}

cd ${system}
obabel ${system}.lig.gast.pdbqt -O ${system}.lig.gast.noH.mol2 -d
mv ${system}.lig.gast.noH.mol2 ${docking_dir}

cd ${docking_dir}
grep '^DOCKED' ${system}.docking.dlg | cut -c9- > ${system}.docking.pdbqt
obabel ${system}.docking.pdbqt -O ${system}.docking.noH.mol2 -d

grep -A8 "LOWEST ENERGY DOCKED CONFORMATION" ${system}.docking.dlg  | tail -n1 | awk '{print $4}' >> index_lowest_score.txt

cd ../../

done
