#!/usr/bin/sh

# This script will take a multimol2 and move the center of mass of each molecule to the origin
# Needs the name of the multimol2 as a command line argument

module load chimera/1.13.1
mol_count="$(grep 'MOLECULE' ${1} | wc -l)"


cat > chimera.com <<EOF
open ${1}
write format mol2  #0 All_Drugs_Chim.mol2
EOF

chimera --nogui chimera.com >& chimera${i}_1.out


rm chimera.com
