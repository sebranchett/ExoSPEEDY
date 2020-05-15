#!/bin/sh

OUT_FILE=test_output.txt
echo 'test results' > $OUT_FILE

# sample the 3 .txt output files
for IN_FILE in ../output/exp_101/atdf101_1979.txt \
               ../output/exp_101/attm101_1979.txt \
               ../output/exp_101/atva101_1979.txt
do
 echo $IN_FILE >> $OUT_FILE

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
IN_FILE=../output/exp_101/atgcm101.lis 
echo $IN_FILE >> $OUT_FILE
tail -8 $IN_FILE >> $OUT_FILE
