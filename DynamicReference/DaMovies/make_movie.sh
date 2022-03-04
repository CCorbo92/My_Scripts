#This script will make growth tree movies given the growth trees were formatted in the way this script expects
#This script assumes the full reference is the first molecule in the multimol2
# You need to change the 3 variables down below

module load anaconda/3
module load chimera/1.13.1

#Start with the growth tree you want to make a movie for in this main directory
file="output.growth_tree_60.mol2"

#Put the pdb code and anchor num here so dont end up with redundant naming of growth tree movies
pdb="PDB_CODE"
anchor="Anchornum"
#This is how many chimera models will have to be loaded for the movie
sizetree=`grep "MOLECULE" $file | wc -l`

mkdir ${file}_dir
cp $file ${file}_dir
cd ${file}_dir

#Run a python script to write a chimera script that is appropriately sized
python3 ../write_chimera_script.py $file $sizetree $pdb $anchor>> chimera_script.com

#Load chimera module and run script to write the movie
