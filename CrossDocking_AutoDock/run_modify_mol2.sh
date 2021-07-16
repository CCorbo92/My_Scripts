system_file="clean.systems.all"

for system in `cat ${system_file}`; do  
cd /gpfs/projects/rizzo/ccorbo/DOCK6_with_ambpdb/SB_2020_testset/${system} 
mkdir zzz.autodockfiles
cd zzz.autodockfiles
python ../../modify_mol2_autodock.py ../${system}.rec.clean.mol2 >> ${system}.rec.clean.mol2
cp ../${system}.lig.am1bcc.mol2 ./${system}.lig.am1bcc.mol2
done
