export system_file=$1
export docking_dir=$2

for system in `  cat ${system_file} `
do
echo ${system}
echo -n  ${system} " " >> /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/PoseReprod_Autodock/21_07_04/timing_all.txt
cd ${system}/${docking_dir}
python /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/PoseReprod_Autodock/21_07_04/timing.py ${system}.docking.dlg >> /gpfs/projects/rizzo/ccorbo/Testing_Grounds/AutoDock/PoseReprod_Autodock/21_07_04/timing_all.txt
cd ../../
done
