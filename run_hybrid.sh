 # g++ hybrid.cc
 g++ `root-config --cflags` hybrid.cc -o hybrid `root-config --libs`
 for Tc in 180 200
  do
  for quench_method in 0 1 2  
   do
   for iRAAfit in 0 1 
    do
    if [ "$iRAAfit" = "0" ]; then
     if [ "$Tc" = "200" ]; then 
      if [ "$quench_method" = "0" ]; then 
       alpha=3.5927
      fi 
      if [ "$quench_method" = "1" ]; then 
       alpha=2.7
      fi
      if [ "$quench_method" = "2" ]; then
       alpha=0.493
      fi
     fi
     if [ "$Tc" = "180" ]; then 
      if [ "$quench_method" = "0" ]; then 
       alpha=2.8225
      fi
      if [ "$quench_method" = "1" ]; then
       alpha=1.9
      fi
      if [ "$quench_method" = "2" ]; then
       alpha=0.42
      fi
     fi
    fi
    if [ "$iRAAfit" = "1" ]; then
     if [ "$Tc" = "200" ]; then
      if [ "$quench_method" = "0" ]; then
       alpha=5.3746
      fi 
      if [ "$quench_method" = "1" ]; then
       alpha=4.35
      fi
      if [ "$quench_method" = "2" ]; then
       alpha=0.6
      fi
     fi
     if [ "$Tc" = "180" ]; then
      if [ "$quench_method" = "0" ]; then
       alpha=4.3232
      fi
      if [ "$quench_method" = "1" ]; then
       alpha=3.1
      fi
      if [ "$quench_method" = "2" ]; then
       alpha=0.51
      fi
     fi
    fi
    # for iCentrality in {0..8}
    for iCentrality in 0
     do
     nohup ./hybrid $quench_method $iCentrality $alpha $Tc 0 > run_{$quench_method}_{$iCentrality}_{$alpha}_{$Tc}_0.out & 
    done
   done
  done
 done