#! /bin/csh
#BSUB -W 240:00
#BSUB -n 12
#BSUB -R span[hosts=1] 
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q gubbins
#BSUB -J P0.3_T0.5
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


set N=2000
set PStar=(0.3)
set TStar=(0.5)
set nbn=(2)
set maxStepSize=(0.1)
set maxVolChange=3

set numSteps=500000000

set grNumSegs=5
set grSegWidth=200.0
set grNumBins=100000
set grBinWidth=0.01
set grInterval=10000000

set rhoNumBins=250000
set rhoBinWidth=0.01
set rhoInterval=10000000

set thermoBlockSize=10000
set printConfigInterv=500000
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
    ./${exe} $N $PStar $TStar ${nbn} ${numSteps} lj ${maxStepSize} ${maxVolChange} ${printConfigInterv} ${thermoBlockSize} ${rhoBinWidth} ${rhoNumBins} ${rhoInterval} ${grSegWidth} ${grNumSegs} ${grBinWidth} ${grNumBins} ${grInterval} `date +'%s'` > out.out
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

