#!/bin/sh 
list_of_fam="smalltest10.txt"
for ref_sys in `cat ${list_of_fam}`; do
echo ${ref_sys} " "${1}
cd ${ref_sys}
cd anc_${1}
/gpfs/projects/rizzo/ccorbo/DOCK_Builds/dock6.21_11_05_Dynamic_Reference_Developmental_V1.6/bin/dock6 -i dn_generic.in  -o layer_passed_dockmol.out
cd ../../
done
