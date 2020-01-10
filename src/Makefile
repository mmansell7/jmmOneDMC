CChen:=mpicxx
CCmac:=/usr/local/bin/gcc-7
CFLAGShen:=-Wall -fopenmp -lgsl -O3
CFLAGSmac:=-Wall -fopenmp -lgsl -O3
DEBUGflags:=-Wall -fopenmp -lgsl -g
EXE:=../bin/jmmOneDMC

default:
	@echo "No target specified."

clean:
	rm -f *.o ${EXE} g*.dat.* rho.dat.* out.out config.dat.* thermo.dat.* errors.* output.*

henrun: henry2
	${EXE} 2000 0.5 0.5 5 2000 lj 0.05 5 100 10 0.01 100000 100 100 5 0.01 100000 100 125 > out.out

henrun2 : henry2
	${EXE} 1000 0.25 0.75 5 10000 lj 0.05 2 10 10 0.01 100000 100 10 5 0.01 100000 100 125 > out.out

henrun3 : henry2
	${EXE} > out.out

debug : readInput.c readInput.h jmmMCState.c jmmMCState.h pot.c pot.h Main.c
	 ${CChen} ${DEBUGflags} -c readInput.c
	 ${CChen} ${DEBUGflags} -c jmmMCState.c
	 ${CChen} ${DEBUGflags} -c pot.c
	 ${CChen} ${DEBUGflags} -c Main.c
	 ${CChen} ${DEBUGflags} Main.o jmmMCState.o readInput.o pot.o -o ${EXE}

henry2 : readInput jmmMCState_henry2 pot Main_henry2
	${CChen} ${CFLAGShen} Main.o jmmMCState.o readInput.o pot.o -o ${EXE}

macbook : jmmMCState.o Main.o pot.o
	${CCmac} ${CFLAGSmac} Main.o jmmMCState.o -o ${EXE}

Main_henry2 : Main.c jmmMCState.o readInput.o pot.o
	${CChen} ${CFLAGShen} -c Main.c

Main_mac : Main.c jmmMCState.h
	${CCmac} ${CFLAGSmac} -c Main.c

jmmMCState_henry2 : jmmMCState.c jmmMCState.h
	${CChen} ${CFLAGShen} -c jmmMCState.c

jmmMCState_macbook : jmmMCState.c jmmMCState.h
	${CCmac} ${CFLAGSmac} -c jmmMCState.c

testInput : readInput
	${CChen} ${CFLAGShen} -c testInput.c
	${CChen} ${CFLAGShen} testInput.o readInput.o jmmMCState.o -o testInput

readInput : jmmMCState_henry2 readInput.c readInput.h
	${CChen} ${CFLAGShen} -c readInput.c

pot : pot.c pot.h
	${CChen} ${CFLAGShen} -c pot.c

compute : compute.cpp compute.h
	${CChen} ${CFLAGShen} -c compute.cpp

.PHONY: run clean
	
