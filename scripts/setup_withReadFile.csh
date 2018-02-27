#! /bin/csh
#BSUB -W 240:00
#BSUB -n 8
#BSUB -R span[hosts=1] 
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q gubbins
#BSUB -J ljM2P4T10
##BSUB -L /bin/csh

set CASESTR=P0.4_T1.0_${LSB_JOBID}
set INFILE=INPUTP4T10

#setenv OMP_NUM_THREADS 1
source /usr/local/apps/mpich3/centos7/intelmpi2016.csh
set headDir=`pwd`/..
set srcDir=${headDir}/src
set JOBDIR=${headDir}/data/lj/${CASESTR}
mkdir -p ${JOBDIR}


##### FOR NON-LSF EXECUTION ############
#set scratchDir=${headDir}/scratch      #
#set LSB_JOBID=123456                   #
########################################


###### FOR LSF EXECUTION ################
#set scratchDir=/scratch/${LSB_JOBID}  #
##set scratchDir=headDir/data/${LSB_JOBID}                                       #
#########################################
#
#mkdir -p ${scratchDir}
#
#cp ./INPUT ${scratchDir}/
#cp ${srcDir}/*.c ${scratchDir}/
#cp ${srcDir}/*.h ${scratchDir}/
#cp ${srcDir}/Makefile ${scratchDir}/
#
#cd ${scratchDir}
#
#echo Compiling ${scratchDir}/Main.c ...
#

cp ./${INFILE} ${JOBDIR}/INPUT
cp ${srcDir}/*.c ${JOBDIR}/
cp ${srcDir}/*.h ${JOBDIR}/
cp ${srcDir}/Makefile ${JOBDIR}/

cd ${JOBDIR}

echo Compiling ${JOBDIR}/Main.c ...



make henry2
echo Compilation completed I think...

./jmmOneDMC > out.out

echo ...Simulation completed I think

#cp ${scratchDir}/* ${JOBDIR}/
#
cd ${headDir}
#rm -r ${scratchDir}

