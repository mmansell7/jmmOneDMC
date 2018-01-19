#! /bin/csh
#BSUB -W 96:00
#BSUB -n 12
#BSUB -R span[hosts=1] 
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q single_chassis
#BSUB -J P0.5_T0.7
##BSUB -L /bin/csh

#setenv OMP_NUM_THREADS 1
source /usr/local/apps/mpich3/centos7/intelmpi2016.csh
set headDir=`pwd`/..


##### FOR NON-LSF EXECUTION ############
#set scratchDir=${headDir}/scratch      #
#set LSB_JOBID=123456                   #
########################################


##### FOR LSF EXECUTION ################
set scratchDir=/scratch/${LSB_JOBID}  #
                                       #
########################################


set srcDir=${headDir}/src
set binDir=${headDir}/Henry2_Debug
mkdir -p ${scratchDir}
cd ${scratchDir}
set OneDNPTMC_srcPath=${srcDir}/Main.c

set exe=OneDNPTMC

set N=2000
set PStar=(0.5)
set TStar=(0.7)

set maxStepSize=(0.1)
set maxVolChange=5

set numSteps=100000000

set grNumSegs=5
set grSegWidth=200.0
set grNumBins=100000
set grBinWidth=0.01
set grInterval=1000000

set rhoNumBins=250000
set rhoBinWidth=0.01
set rhoInterval=1000000

set thermoBlockSize=10000
set printConfigInterv=100000
#echo `pwd`
cp ${OneDNPTMC_srcPath} ${scratchDir}/Main.c
#echo Compiling a file from ${srcDir}
echo Compiling ${scratchDir}/Main.c ...
mpicc -Wall -fopenmp -lgsl -O3 Main.c -o ./${exe}
echo Compilation completed I think...
#cp ${binDir}/${exe} ${scratchDir}/

#foreach ii ( 1 )
#for ii in 0 1 2 3 4 5 6 7 8 9 10 11 12; do
#for lStarBar in 1.0 1.5 2.0 2.5 3.0 4.0 5.0 7.0 10.0 15.0 20.0 25.0; do
  #echo $lStarBar;
  #maxStepSize=$(bc <<< "scale=4;$lStarBar/2")
  #foreach TStar ( 0.3 )
  #for TStar in 0.1 0.2 0.3 0.5 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10.0 15.0 20.0; do
    #echo Submitting simulation for: $N ${lStarBar[${ii}]} ${TStar} ${numSteps} ${thermoBlockSize} ${maxStepSize[${ii}]} ${rcutStar} ${grNumBins} ${grBlockSize} ${printConfigInterv}
    echo Submitting simulation...
    #./${LJexe} $N ${lStarBar[${ii}]} ${TStar} ${numSteps} ${thermoBlockSize} ${maxStepSize[${ii}]} ${rcutStar} ${grNumBins} ${grBlockSize} ${printConfigInterv} > l${lStarBar[${ii}]}_T${TStar}_${LSB_JOBID}.out
    ./${exe} $N $PStar $TStar ${numSteps} lj ${maxStepSize} ${maxVolChange} ${printConfigInterv} ${thermoBlockSize} ${rhoBinWidth} ${rhoNumBins} ${rhoInterval} ${grSegWidth} ${grNumSegs} ${grBinWidth} ${grNumBins} ${grInterval} > out.out
    #./${exe}
    echo ...Simulation completed I think

    set JOBDIR=${headDir}/data/P${PStar}_T${TStar}_${LSB_JOBID}
    mkdir -p ${JOBDIR}
    cp ${scratchDir}/* ${JOBDIR}/
    #rm ${scratchDir}/*.out ${scratchDir}/gr.dat
    
#  end
#end

cd ${headDir}
rm -r ${scratchDir}

