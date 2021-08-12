sort -n  PDB_ProteinName.csv >  PDB_ProteinName_sort.csv
paste -d ','  Descriptors_sort.csv  FDA_Descriptors.txt >>temp1.csv
paste -d ',' temp1.csv PDB_ProteinName_sort.csv >> temp2.csv
paste -d ',' temp2.csv CrossDock_All_sys.csv >> Full_spreadsheet.csv
rm temp*
