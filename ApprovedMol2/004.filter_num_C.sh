mkdir Removed
cd RDKit
grep "CCCCCCCC" * | awk -F':' '{print $1}' >> Failure_list.txt
while read line; do
mv ${line} ../Removed
done<"Failure_list.txt"
rm Failure_list.txt
