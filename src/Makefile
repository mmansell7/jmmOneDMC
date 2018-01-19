CChen:=mpicc
CCmac:=/usr/local/bin/gcc-7
CFLAGShen:=-Wall -fopenmp -lgsl -O3
CFLAGSmac:=-Wall -fopenmp -lgsl -O3
EXE:=jmmOneDMC

default:
	@echo "No target specified."

clean:
	rm -f *.o ${EXE}

run: ${EXE}
	./${EXE}

henry2 : jmmMCState.o Main.o
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
	
