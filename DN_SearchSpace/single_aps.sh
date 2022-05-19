for file in `cat list.txt`; do
   dummy_count=`grep -wc "Du" $file`
   if [ $dummy_count -gt 1 ]
   then
      start=2
      count=$(( 2*${dummy_count} ))
      for (( x=$start; x<=$count; x++ )); do
         a=2
         b=4
         c=6
         d=8

         if [ $x == 2 ]
	    then
		a=10
	 fi

         if [ $x == 4 ]
            then
                b=10
         fi

         if [ $x == 6 ]
            then
                c=10
         fi

         if [ $x == 8 ]
            then
                d=10
         fi

         echo $x
         echo $a " " $b " " $c " " $d " " $e

         sed ":a;N;\$!ba;s/Du/H /$d; s/Du/H /$c; s/Du/H /$b; s/Du/H /$a" $file >> ./processed/${x}_${file}

         x=$(( ${x}+1 ))
      done
   fi
done
