# [method][Tc][RAAfit]

 g++ hybrid.cc
 
 for Tc in 180 200
  do
  for quench_method in 1 4 3 
   do
   for iRAAfit in 0 1 
    do
    if [ "$iRAAfit" = "0" ]; then
     if [ "$Tc" = "200" ]; then 
      if [ "$quench_method" = "1" ]; then 
       alpha=3.5927
      fi 
      if [ "$quench_method" = "4" ]; then 
       alpha=2.629
      fi
      if [ "$quench_method" = "3" ]; then
       alpha=0.4759
      fi
     fi
     if [ "$Tc" = "180" ]; then 
      if [ "$quench_method" = "1" ]; then 
       alpha=5.3746
      fi
      if [ "$quench_method" = "4" ]; then
       alpha=4.3232
      fi
      if [ "$quench_method" = "3" ]; then
       alpha=4.1035
      fi
     fi
    fi
    if [ "$iRAAfit" = "0" ]; then
     if [ "$Tc" = "200" ]; then
      if [ "$quench_method" = "1" ]; then
       alpha=5.3746
      fi 
      if [ "$quench_method" = "4" ]; then
       alpha=4.1035
      fi
      if [ "$quench_method" = "3" ]; then
       alpha=0.5709
      fi
     fi
     if [ "$Tc" = "180" ]; then
      if [ "$quench_method" = "1" ]; then
       alpha=4.3232
      fi
      if [ "$quench_method" = "4" ]; then
       alpha=2.8127
      fi
      if [ "$quench_method" = "3" ]; then
       alpha=0.4714
      fi
     fi
    fi
    # for iCentrality in {0..8}
    for iCentrality in 0
     do
     ./a.out $quench_method $iCentrality $alpha $Tc
    done
   done
  done
 done