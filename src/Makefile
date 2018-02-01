CChen:=mpicc
CCmac:=/usr/local/bin/gcc-7
CFLAGShen:=-Wall -fopenmp -lgsl -O3
CFLAGSmac:=-Wall -fopenmp -lgsl -O3
EXE:=jmmOneDMC

default:
	@echo "No target specified."

clean:
	rm -f *.o ${EXE} g*.dat rho.dat out.out config.dat thermo.dat

henrun: henry2
	./${EXE} 2000 0.5 0.5 5 2000 lj 0.05 5 100 10 0.01 100000 100 100 5 0.01 100000 100 125 > out.out

henrun2 : henry2
	./${EXE} 1000 0.25 0.75 5 10000 lj 0.05 2 10 10 0.01 100000 100 10 5 0.01 100000 100 125 > out.out

henry2 : Main_henry2 jmmMCState_henry2 jmmMCState.o Main.o
	${CChen} ${CFLAGShen} Main.o jmmMCState.o -o ./${EXE}

macbook : jmmMCState.o Main.o
	${CCmac} ${CFLAGSmac} Main.o jmmMCState.o -o ./${EXE}

Main_henry2 : Main.c jmmMCState.h
	${CChen} ${CFLAGShen} -c Main.c

Main_mac : Main.c jmmMCState.h
	${CCmac} ${CFLAGSmac} -c Main.c

jmmMCState_henry2 : jmmMCState.c jmmMCState.h
	${CChen} ${CFLAGShen} -c jmmMCState.c

jmmMCState_macbook : jmmMCState.c jmmMCState.h
	${CCmac} ${CFLAGSmac} -c jmmMCState.c

.PHONY: run clean
	
