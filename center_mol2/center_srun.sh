#!/usr/bin/sh

# This script will take a multimol2 and move the center of mass of each molecule to the origin
# Needs the name of the multimol2 as a command line argument


#################################
# Separate into individual mol2
cat > chimera${2}.com <<EOF
open ${1}
write format mol2  #0.${2} ${2}.mol2
EOF

chimera --nogui chimera${2}.com >& chimera${2}_1.out

rm chimera${2}.com
################################
# Find the x,y,z center of mass for each mol
cat > chimera${2}.com <<EOF
open ${2}.mol2
define centroid mass true radius 0.0
EOF

chimera --nogui chimera${2}.com >& chimera${2}_2.out

#################################
# If the line in the output has 10 columns
num_col="$(grep "centroid" chimera${2}_2.out | awk '{print NF}')"

if [ "${num_col}" == "10" ]; then
  x="$(grep "centroid name" chimera${2}_2.out | awk '{print $8}' | tr -d ',')"
  y="$(grep "centroid name" chimera${2}_2.out | awk '{print $9}' | tr -d ',')"
  z="$(grep "centroid name" chimera${2}_2.out | awk '{print $10}' | tr -d ')')"

    if  echo ${x} | grep -q "-" 
      then
        x_temp="$(echo ${x} | tr -d '-')"
        x="${x_temp}"
    else
        x_temp="-${x}"
        x="${x_temp}"
    fi

    if  echo ${y} | grep -q "-" 
      then
        y_temp="$(echo ${y} | tr -d '-')"
        y="${y_temp}"
    else
        y_temp="-${y}"
        y="${y_temp}"
    fi

    if  echo ${z} | grep -q "-" 
      then
        z_temp="$(echo ${z} | tr -d '-')"
        z="${z_temp}"
    else
        z_temp="-${z}"
        z="${z_temp}"
    fi

fi ############################# 10 columns

#################################
# If the line in the output has 9 columns
if [ "${num_col}" == "9" ]; then
x="$(grep "centroid name" chimera${2}_2.out | awk '{print $7}' | tr -d ','| tr -d '(')"
y="$(grep "centroid name" chimera${2}_2.out | awk '{print $8}' | tr -d ',')"
z="$(grep "centroid name" chimera${2}_2.out | awk '{print $9}' | tr -d ')')"

  if  echo ${x} | grep -q "-"
    then
        x_temp="$(echo ${x} | tr -d '-')"
        x="${x_temp}"
  else
        x_temp="-${x}"
        x="${x_temp}"
  fi

  if  echo ${y} | grep -q "-"
    then
        y_temp="$(echo ${y} | tr -d '-')"
        y="${y_temp}"
  else
        y_temp="-${y}"
        y="${y_temp}"
  fi

  if  echo ${z} | grep -q "-"
    then
        z_temp="$(echo ${z} | tr -d '-')"
        z="${z_temp}"
  else
        z_temp="-${z}"
        z="${z_temp}"
  fi

fi ############################# 9 columns
rm chimera${2}.com

################################
# translate the molecule to origin
cat > chimera${2}.com <<EOF
open ${2}.mol2
move ${x},${y},${z}
write format mol2  #0 ${2}_centered.mol2
EOF

chimera --nogui chimera${2}.com >& chimera${2}_3.out

rm chimera${2}.com

mv chimera${2}*out  ./outfiles
rm ${2}.mol2

