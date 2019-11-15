#!/bin/ksh

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file ( 0 = no restart ) 

# Define directory names
# set -x

if [ $# -ne 3 ] ; then

    echo 'Usage: '$0' resol. exp_no. restart_no' 1>&2
    exit 1

fi

UT=..	
SA=$UT/source
#mkdir $UT/tmp_$2	
#CA=$UT/tmp_$2
CA=$UT/tmp
mkdir $UT/output/exp_$2	
CB=$UT/output/exp_$2
CC=$UT/ver41.input
CD=$UT/output/exp_$3	

# Edit input files if needed and store them with exp code

echo "Do you want to modify the time-stepping parameters (y/n)?"
read MODIFY
if [ $MODIFY = 'y' ] ; then
  emacs $UT/ver41.5.input/cls_instep.h $SA/doc_instep.txt
fi

echo "Do you want to modify the dynamics parameters (y/n)?"
read MODIFY
if [ $MODIFY = 'y' ] ; then
  emacs $UT/ver41.5.input/cls_indyns.h  $SA/doc_indyns.txt
fi

echo "Do you want to modify the physics parameters (y/n)?"
read MODIFY
if [ $MODIFY = 'y' ] ; then
  emacs $UT/ver41.5.input/cls_inphys.h  $SA/doc_inphys.txt
fi

echo "Do you want to modify the land model parameters (y/n)?"
read MODIFY
if [ $MODIFY = 'y' ] ; then
  emacs $UT/ver41.5.input/cls_inland.h 
fi

echo "Do you want to modify the sea/ice model parameters (y/n)?"
read MODIFY
if [ $MODIFY = 'y' ] ; then
  emacs $UT/ver41.5.input/cls_insea.h
fi
	
echo "Do you want to modify the input files (y/n)?"
read MODIFY
if [ $MODIFY = 'y' ] ; then
  emacs $UT/ver41.5.input/inpfiles.s
fi


mkdir $UT/input/exp_$2

echo "model version   :   41"  > $UT/input/exp_$2/run_setup
echo "hor. resolution : " $1  >> $UT/input/exp_$2/run_setup
echo "experiment no.  : " $2  >> $UT/input/exp_$2/run_setup
echo "restart exp. no.: " $3  >> $UT/input/exp_$2/run_setup
	
# Copy files from basic version directory

echo "copying from $SA/source to $CA"
rm -f $CA/*

cp $SA/makefile $CA/
cp $SA/*.f      $CA/
cp $SA/*.h      $CA/
cp $SA/*.s      $CA/

cp $CA/par_horres_$1.h   $CA/atparam.h
cp $CA/par_verres.h      $CA/atparam1.h 

# Copy parameter and namelist files from user's .input directory

echo "ver41.input new files ..."
ls $UT/ver41.input

echo "copying parameter and namelist files from $UT/ver41.input "
cp $UT/ver41.5.input/cls_*.h     $CA/
cp $UT/ver41.5.input/inpfiles.s  $CA/
cp $UT/ver41.5.input/cls_*.h     $UT/input/exp_$2
cp $UT/ver41.5.input/inpfiles.s  $UT/input/exp_$2

# Copy modified model files from user's update directory

echo "update new files ..."
ls $UT/update

echo "copying modified model files from $UT/update"
cp $UT/update/*.f   $CA/
cp $UT/update/*.h   $CA/
cp $UT/update/make* $CA/	
cp $UT/update/*.f   $UT/input/exp_$2
cp $UT/update/*.h   $UT/input/exp_$2
cp $UT/update/make* $UT/input/exp_$2
			
# Set input files

cd $CA

# Set experiment no. and restart file (if needed)

echo $3 >  fort.2
echo $2 >> fort.2

if [ $3 != 0 ] ; then
  echo "link restart file atgcm$3.rst to fort.3"
  ln -s $CD/atgcm$3.rst fort.3
fi 

# Link input files

echo 'link input files to fortran units'

ksh inpfiles.s $1

ls -l fort.*

echo ' compiling at_gcm - calling make'

make imp.exe  

#
# create and execute a batch job to run the model
#

cat > run.job << EOF1
set -x
 
cd $CA
pwd


echo 'the executable file...'
ls -l imp.exe

 
time ./imp.exe > out.lis

mv out.lis $CB/atgcm$2.lis
mv fort.10 $CB/atgcm$2.rst

mv at*$2.ctl   $CB
mv at*$2_*.grd $CB

mv day*$2.ctl   $CB
mv day*$2_*.grd $CB

cd $CB

chmod 644 at*$2.* 

EOF1


ksh run.job &

exit
