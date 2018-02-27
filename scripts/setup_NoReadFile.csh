#! /bin/csh
#BSUB -W 240:00
#BSUB -n 12
#BSUB -R span[hosts=1] 
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q gubbins
#BSUB -J P1.0_T0.5
##BSUB -L /bin/csh

#setenv OMP_NUM_THREADS 1
source /usr/local/apps/mpich3/centos7/intelmpi2016.csh
set headDir=`pwd`/..
set srcDir=${headDir}/src
set OneDNPTMC_srcPath=${srcDir}/Main.c
set exe=jmmOneDMC


##### FOR NON-LSF EXECUTION ############
#set scratchDir=${headDir}/scratch      #
#set LSB_JOBID=123456                   #
########################################


##### FOR LSF EXECUTION ################
set scratchDir=/scratch/${LSB_JOBID}  #
                                       #
########################################

mkdir -p ${scratchDir}
cd ${scratchDir}

set srcDir=${headDir}/src
#set OneDNPTMC_srcPath=${srcDir}/Main.c
set binDir=${headDir}/bin
set exe=jmmOneDMC

set N=2000
set PStar=(1.0)
set TStar=(0.5)
set nbn=(2)
set maxStepSize=(0.1)
set maxVolChange=3

set numSteps=5000000000

set grNumSegs=10
set grSegWidth=200.0
set grNumBins=100000
set grBinWidth=0.01
set grInterval=10000000

set rhoNumBins=250000
set rhoBinWidth=0.01
set rhoInterval=10000000

set thermoBlockSize=10000

set printConfigInterv=500000

cp ${OneDNPTMC_srcPath} ${scratchDir}/Main.c
cp ${srcDir}/jmmMCState.c ${scratchDir}/
cp ${srcDir}/jmmMCState.h ${scratchDir}/
cp ${srcDir}/Makefile ${scratchDir}/
pwd
ls -hl

echo Compiling ${scratchDir}/Main.c ...
make henry2
#mpicc -Wall -fopenmp -lgsl -O3 Main.c -o ./${exe}
echo Compilation completed I think...

./${exe} $N $PStar $TStar ${nbn} ${numSteps} lj ${maxStepSize} ${maxVolChange} ${printConfigInterv} ${thermoBlockSize} ${rhoBinWidth} ${rhoNumBins} ${rhoInterval} ${grSegWidth} ${grNumSegs} ${grBinWidth} ${grNumBins} ${grInterval} `date +'%s'` > out.out

echo ...Simulation completed I think

set JOBDIR=${headDir}/data/P${PStar}_T${TStar}_${LSB_JOBID}
mkdir -p ${JOBDIR}
cp ${scratchDir}/* ${JOBDIR}/

cd ${headDir}
rm -r ${scratchDir}

