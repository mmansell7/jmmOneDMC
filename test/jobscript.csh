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

set EXE=/gpfs_backup/gubbins_data/jmmansel/Applications/jmmOneDMC/src/jmmOneDMC

rm -f config.dat.mcs thermo.dat.mcs rho.dat.mcs g0.dat.mcs g1.dat.mcs g2.dat.mcs g3.dat.mcs g4.dat.mcs g5.dat.mcs g6.dat.mcs g7.dat.mcs g8.dat.mcs g9.dat.mcs out.out

echo "Executing program" ${EXE}
${EXE} > out.out
echo "Program" ${EXE} "completed"

