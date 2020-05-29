#!/bin/sh
set -e

# check a sample of values in .txt output files
diff -w checked_ocean_output.txt test_ocean_output.txt

# check final values of .lis file
prec=".0001"
A="../output/exp_200/atgcm200.lis"
B="../output/exp_201/atgcm201.lis"

for ROW in 374 380
do
  for COLUMN in {4..6}
  do
    checked=$(awk -v row="$ROW" -v column="$COLUMN" 'NR == row {print $column}' ${A})
    test=$(awk -v row="$ROW" -v column="$COLUMN" 'NR == row {print $column}' ${B})
    
    # check if difference is larger than 1%
    toobig=$(echo "(1. - $test / $checked)^2 > $prec" | bc -l)
    if [ $toobig -ne 0 ] ; then echo "error too large line $ROW, field $COLUMN"; exit 1; fi 
  done
done
