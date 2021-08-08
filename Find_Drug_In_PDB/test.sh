while read name;
do
  wget https://drugcentral.org/?q=${name}
  echo -n "drugcentral.org" >> temp_url.txt
  grep -m1 "drugcard" index.html\?q\=${name} | cut -c 10-| rev | cut -c 3- | rev >> temp_url.txt
  while read line;
  do
  wget ${line}
  done < "temp_url.txt"
  cat temp_url.txt | cut -c 26- >> temp_url2.txt
  while read line;
  do
  echo -n ${name} >> target_class.txt
  echo -n " " >> target_class.txt
  grep -A20 "Target</th><th>Class" ${line} | grep "GPCR" | head -n1 | tr -d "td" >> target_class.txt
  grep -A20 "Target</th><th>Class" ${line} | grep "Ion channel" | head -n1 | tr -d "td" >> target_class.txt
  grep -A20 "Target</th><th>Class" ${line} | grep "Kinase"| head -n1 |  tr -d "td" >> target_class.txt
  grep -A20 "Target</th><th>Class" ${line} | grep "Nuclear hormone receptor" | head -n1 | tr -d "td" >> target_class.txt
  echo " " >> target_class.txt
  done < "temp_url2.txt"
  #rm temp*
  #rm *q=*
done < "test.txt"
