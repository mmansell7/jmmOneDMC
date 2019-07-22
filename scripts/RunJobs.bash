#!/bin/bash

# Uncomment lines below to run with LJ potential
pot=LJ
potabr=LJ
m=-1

# Uncomment lines below to run with harmonic potential
#pot=HARMONIC
#potabr=HARM
#m=-1

curDir=`pwd`

N=80
Pressures=( 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 )
Temperatures=( 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 )

for P in "${Pressures[@]}"
do
   for T in "${Temperatures[@]}"
   do

      echo Pot: ${pot}  m: ${m}   P: ${P}   T: ${T}

      f=../data/${pot}/m${m}/N${N}/P${P}_T${T}
      if ls ${f}* 1> /dev/null 2>&1; then
         echo Directory exists.  Skipping.
      else
         echo Directory does not exist.  Preparing this case.
         inFile=${f}/INPUT
         jobFile=${f}/job.csh
         echo Making directory ${f}
         mkdir -p ${f}
         echo Making input file ${inFile}


cat >${inFile} << EOL
N          ${N}
P          ${P}
T          ${T}
RELAX
NUMSTEPS   5000000000
POT        ${pot}
NBN        ${m}
MAXSTEP    0.1
MAXDV      2.0
CPI        50000000
TPI        10000000
RBW        0.1
RHONB      1000
RHOPI      50000000
GSW        200.0
GNS        10
GBW        0.1
GNB        1000
GPI        50000000
SEED       92847
ENGCHECK   10000
DADJ       1000000
VADJ       1000000
EOL


         echo Making job script file ${jobFile}


cat >${jobFile} << EOL
#! /bin/csh
#BSUB -W 80:00
#BSUB -n 8
#BSUB -R span[hosts=1] 
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q single_chassis
#BSUB -J ${potabr}M${m}N${N}P${P}T${T}
##BSUB -L /bin/csh

set CASESTR=P${P}_T${T}_m${m}_\${LSB_JOBID}

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
         cd ${f}
         bsub < job.csh
         cd ${curDir}

      fi
   done
done




