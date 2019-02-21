#!/bin/bash

# Uncommented lines below to run with LJ potential
#pot=LJ
#potabr=LJ
#m=-1

# Uncomment lines below to run with harmonic potential
pot=HARMONIC
potabr=HARM
m=-1

curDir=`pwd`

Temperatures=( 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.2 2.4 2.6 2.8 3.0 5.0 10.0 )
#Pressures=( 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 5.0 10.0 )
Pressures=( 5.0 10.0 )

for P in "${Pressures[@]}"
do
   for T in "${Temperatures[@]}"
   do
#Pressures=(    0.6 )
#Temperatures=( 0.5 )
#for ii in `seq 0 41`
#do
#      P=${Pressures[$ii]}
#      T=${Temperatures[$ii]}

      echo Pot: ${pot}  m: ${m}   P: ${P}   T: ${T}

      f=../data/${pot}/m${m}/Running_Jobs/P${P}_T${T}
      if ls ../data/${pot}/m${m}/Running_Jobs/P${P}_T${T}_* 1> /dev/null 2>&1; then
         echo Directory exists.  Skipping.
      else
         echo Directory does not exist.  Preparing this case.
         newDir=../data/${pot}/m${m}/Running_Jobs/P${P}_T${T}_m${m}
         echo Making directory ${newDir}
         mkdir -p ${newDir}
         inFile=${newDir}/INPUT
         echo Making input file ${inFile}


cat >${inFile} << EOL
N          2000
P          ${P}
T          ${T}
RELAX
NUMSTEPS   5000000000
POT        ${pot}
NBN        ${m}
MAXSTEP    0.1
MAXDV      2.0
CPI        500000
TPI        10000
RBW        0.1
RHONB      1000
RHOPI      1000000
GSW        200.0
GNS        10
GBW        0.1
GNB        1000
GPI        1000000
SEED       92847
ENGCHECK   10000
DADJ       100000
VADJ       100000
EOL


         jobFile=${newDir}/job.csh
         echo Making job script file ${jobFile}


cat >${jobFile} << EOL
#! /bin/csh
#BSUB -W 80:00
#BSUB -n 8
#BSUB -R span[hosts=1] 
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q single_chassis
#BSUB -J ${potabr}M${m}P${P}T${T}
##BSUB -L /bin/csh

set CASESTR=P${P}_T${T}_m${m}\${LSB_JOBID}

source /usr/local/apps/mpich3/centos7/intelmpi2016.csh
set srcDir=/gpfs_backup/gubbins_data/jmmansel/jmmOneDMC/src/
cp \${srcDir}/*.c ./
cp \${srcDir}/*.h ./
cp \${srcDir}/Makefile ./

echo Compiling Main.c ...
make henry2
echo Compilation completed I think...

./jmmOneDMC > out.out

echo ...Simulation completed I think
EOL


         echo Submitting job.
         cd ${newDir}
         bsub < job.csh
         cd ${curDir}

      fi
   done
done




