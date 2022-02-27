#This script will set up viewdock states automatically with dynamic reference molecules and standard molecules compared to full reference

file="../smalltest20.txt"
directory="2022_02_22_20SYS_80Anchors_RMSDCluster_0.2"

cd $directory

for pdb in `cat $file`;
do

cat > chimera_${pdb}.com << EOF
viewdock C:/Users/ccxaf/Downloads/RizzoLab/DynamicRef/PlottingData/${directory}/All_Results/${pdb}_rescored_fullref_25_ranked.mol2
viewdock C:/Users/ccxaf/Downloads/RizzoLab/DynamicRef/PlottingData/${directory}/All_Results_NonDynamic/${pdb}_rescored_fullref_25_ranked.mol2
open C:/Users/ccxaf/Downloads/RizzoLab/DynamicRef/ReferenceSystems/${pdb}.lig.am1bcc.mol2

background solid white

delete HC

color dim gray #2
color byhet #2
setattr m stickScale 0.45 #2

color orange #0
color byhet #0

color green #1
color byhet #1

save C:/Users/ccxaf/Downloads/RizzoLab/DynamicRef/PlottingData/${directory}/${pdb}.py

EOF

done
