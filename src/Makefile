GIT_VERSION := "$(shell git describe --abbrev=40 --always)"

CChen:=mpicxx
CCmac:=/usr/local/bin/gcc-7
CFLAGShen:=-Wall -fopenmp -lgsl -O3 -D'VERSION=$(GIT_VERSION)'
CFLAGSmac:=-Wall -fopenmp -lgsl -O3 -D'VERSION=$(GIT_VERSION)'
DEBUGflags:=-Wall -fopenmp -lgsl -g -D'VERSION=$(GIT_VERSION)'
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

debug : readInput.cpp readInput.h jmmMCState.cpp jmmMCState.h pot.cpp pot.h Main.cpp compute.h compute.cpp
	${CChen} ${DEBUGflags} -c readInput.cpp
	${CChen} ${DEBUGflags} -c jmmMCState.cpp
	${CChen} ${DEBUGflags} -c pot.cpp
	${CChen} ${DEBUGflags} -c Main.cpp
	${CChen} ${DEBUGflags} -c compute.cpp
	${CChen} ${DEBUGflags} Main.o jmmMCState.o readInput.o pot.o compute.cpp -o ${EXE}

henry2 : readInput jmmMCState_henry2 pot Main_henry2 compute
	${CChen} ${CFLAGShen} Main.o jmmMCState.o readInput.o pot.o compute.o -o ${EXE}

macbook : jmmMCState.o Main.o pot.o
	${CCmac} ${CFLAGSmac} Main.o jmmMCState.o -o ${EXE}

Main_henry2 : Main.cpp jmmMCState.o readInput.o pot.o
	${CChen} ${CFLAGShen} -c Main.cpp

Main_mac : Main.cpp jmmMCState.h
	${CCmac} ${CFLAGSmac} -c Main.cpp

jmmMCState_henry2 : jmmMCState.cpp jmmMCState.h
	${CChen} ${CFLAGShen} -c jmmMCState.cpp

jmmMCState_macbook : jmmMCState.cpp jmmMCState.h
	${CCmac} ${CFLAGSmac} -c jmmMCState.cpp

testInput : readInput
	${CChen} ${CFLAGShen} -c testInput.c
	${CChen} ${CFLAGShen} testInput.o readInput.o jmmMCState.o -o testInput

readInput : jmmMCState_henry2 readInput.cpp readInput.h
	${CChen} ${CFLAGShen} -c readInput.cpp

pot : pot.cpp pot.h
	${CChen} ${CFLAGShen} -c pot.cpp

compute : compute.cpp compute.h
	${CChen} ${CFLAGShen} -c compute.cpp

.PHONY: run clean
	
