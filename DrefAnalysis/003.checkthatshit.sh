grep "NO ANCHOR MATCHES" */*/layer_passed_dockmol.out  | awk -F'/' '{print $1 " " $2}' > Missed_Anchors.txt
grep "It has a score of" */*/layer_passed_dockmol.out | awk '{print $6}' > Best_Anchor_Match_HMS.txt
grep "It has a score of" */*/layer_passed_dockmol.out | awk -F'/' '{print $1 " " $2}' > Definitely_Successful.txt
grep "The root size after anchor orienting is 0" */*/layer_passed_dockmol.out | awk -F'/' '{print $1 " " $2}' > Orientation_Failure.txt
#cat Missed_Anchors.txt Orientation_Failure.txt Definitely_Successful.txt | sort -n > AccountedFor.txt
#awk '{print $1}'  AccountedFor.txt  | uniq > Systemshere.txt
#diff sortFARMA.txt Systemshere.txt | grep "<" > Missing_Systems.txt
