while read line num;
do 
grep ${line} Spreadsheet.csv | awk '{print $1}' >> ${line}.txt
done <  "zzz.CrossDockingFamilies.txt"
