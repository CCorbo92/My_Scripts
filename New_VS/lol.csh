
set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set mpidir    = "${VS_MPIDIR}/bin"

set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set descriptdir = "${rootdir}/zzz.descriptor"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"
cd ${rootdir}/${system}/007.cartesian-min-and-rescore/${vendor}/
set num_chunks = "9"
set chunk = "0" 
mkdir chunks
while (${chunk} < ${num_chunks})
cd ${rootdir}/${system}/007.cartesian-min-and-rescore/${vendor}/
mv chunk_${chunk}/chunk_${chunk}_cart_min_rescore_scored.mol2 chunks
@ chunk++
end 
