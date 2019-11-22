#! /bin/csh
#BSUB -W 2:00
##BSUB -n 8
##BSUB -R span[hosts=1] 
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q single_chassis
#BSUB -J jmmOneDMCTest
##BSUB -L /bin/csh

source /usr/local/apps/mpich3/centos7/intelmpi2016.csh
set EXE=/gpfs_backup/gubbins_data/jmmansel/jmmOneDMC/src/jmmOneDMC

rm -f config.dat.mcs thermo.dat.mcs rho.dat.mcs g*.dat.mcs

echo "Executing program" ${EXE}
${EXE} > out.out
echo "Program" ${EXE} "completed"

