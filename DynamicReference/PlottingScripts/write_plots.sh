#One main directory for experiment and two subdirectory for dynamic and standard
cd 2022_02_14_80Anchors_DrefStd

echo "Merge all values for entire distribution"
for desc in `cat ../Descriptors_List.txt`;
do
        echo "25"
	cat All_Results/*_${desc}_25.txt > All_Results/ALL_${desc}_25.txt
	cat All_Results_NonDynamic/*_${desc}_25.txt > All_Results_NonDynamic/ALL_${desc}_25.txt

	echo "100" 
	cat All_Results/*_${desc}_100.txt > All_Results/ALL_${desc}_100.txt
        cat All_Results_NonDynamic/*_${desc}_100.txt > All_Results_NonDynamic/ALL_${desc}_100.txt
	
	echo "ALL"
        cat All_Results/*_${desc}_ALL.txt > All_Results/ALL_${desc}_ALL.txt
	cat All_Results_NonDynamic/*_${desc}_ALL.txt > All_Results_NonDynamic/ALL_${desc}_ALL.txt

done
#List of PDB Codes here PLot descriptors for each pdb sys
echo "Generating plots"
for sys in `cat ../smalltest20.txt`;
do 
        echo $sys " ALL"
	python3 ../DN_Plotting_Property_Distribution.py $sys  ../Descriptors_List.txt ALL
done

for sys in `cat ../smalltest20.txt`;
do
        echo $sys " 100"
        python3 ../DN_Plotting_Property_Distribution.py $sys  ../Descriptors_List.txt 100
done

for sys in `cat ../smalltest20.txt`;
do
        echo $sys " 25"
        python3 ../DN_Plotting_Property_Distribution.py $sys  ../Descriptors_List.txt 25
done
