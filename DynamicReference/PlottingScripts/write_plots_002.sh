#One main directory for experiment and two subdirectory for dynamic and standard
cd 2022_02_22_20SYS_80Anchors_RMSDCluster_0.2

echo "Merge all values for entire distribution"
rm All_Results*/ALL*ALL.txt
for desc in `cat ../Descriptors_List.txt`;
do
	if [ "$desc" = "GRD" ] || [ "$desc" = "HMS" ]; then
		head -n 25 All_Results/*${desc}_ALL.txt | grep -v ">" | grep -v -e '^$' >  All_Results/${desc}_ALL_25.txt
        	head -n 25 All_Results_NonDynamic/*${desc}_ALL.txt | grep -v ">" | grep -v -e '^$'  >  All_Results_NonDynamic/${desc}_ALL_25.txt	

		head -n 100 All_Results/*${desc}_ALL.txt | grep -v ">" | grep -v -e '^$' >  All_Results/${desc}_ALL_100.txt
		head -n 100 All_Results_NonDynamic/*${desc}_ALL.txt | grep -v ">" | grep -v -e '^$'  >  All_Results_NonDynamic/${desc}_ALL_100.txt
	fi 
        if [ "$desc" = "VOL" ] || [ "$desc" = "TAN" ]; then
		tail -n 25 All_Results/*${desc}_ALL.txt | grep -v ">" | grep -v -e '^$' >  All_Results/${desc}_ALL_25.txt
		tail -n 25 All_Results_NonDynamic/*${desc}_ALL.txt | grep -v ">" | grep -v -e '^$'  >  All_Results_NonDynamic/${desc}_ALL_25.txt

                tail -n 100 All_Results/*${desc}_ALL.txt | grep -v ">" | grep -v -e '^$' >  All_Results/${desc}_ALL_100.txt
		tail -n 100 All_Results_NonDynamic/*${desc}_ALL.txt | grep -v ">" | grep -v -e '^$'  >  All_Results_NonDynamic/${desc}_ALL_100.txt
	fi

done
#List of PDB Codes here PLot descriptors for each pdb sys
echo "Generating plots"
python3 ../DN_Plotting_Property_Distribution_with_specified_number.py   ../Descriptors_List.txt 25
python3 ../DN_Plotting_Property_Distribution_with_specified_number.py   ../Descriptors_List.txt 100
