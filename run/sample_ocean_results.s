#!/bin/sh

OUT_FILE=test_ocean_output.txt
echo 'test results' > $OUT_FILE

# sample the 3 .txt output files
for IN_FILE in ../output/exp_201/atdf201_2000.txt \
               ../output/exp_201/attm201_2000.txt \
               ../output/exp_201/atva201_2000.txt
do
 echo '--------' >> $OUT_FILE

 # find length of the results file and the middle
 LENGTH=$(wc -l $IN_FILE | awk '{print $1}')
 let "MIDDLE = $LENGTH / 2"
 wc -l $IN_FILE | awk '{print $1}' >> $OUT_FILE

 # capture the start, the middle and the end
 head -5 $IN_FILE >> $OUT_FILE
 head -$MIDDLE $IN_FILE | tail -5 >> $OUT_FILE
 tail -5 $IN_FILE >> $OUT_FILE
done

# capture the final results in the .lis output file
IN_FILE=../output/exp_201/atgcm201.lis 
echo '--------' >> $OUT_FILE
tail -8 $IN_FILE >> $OUT_FILE
