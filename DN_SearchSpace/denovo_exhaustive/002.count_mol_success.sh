cd runs
wc -l */count.txt | awk '{print $1}' | grep -w  -v "0" | wc -l
wc -l */count.txt | awk '{print $1}' | grep -w  -v "0"
