#!/bin/sh 
ref_sys="1MES"
cd ${ref_sys}
cd anc_${1}
/gpfs/projects/rizzo/ccorbo/dock6.1_16_21_double_end_postprocess_new_beta/bin/dock6 -i dn_generic.in -o dn_generic.out

