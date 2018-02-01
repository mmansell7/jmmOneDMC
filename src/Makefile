CChen=mpicc
CCmac=/usr/local/bin/gcc-7
CFLAGShen=-Wall -fopenmp -lgsl -O3
CFLAGSmac=-Wall -fopenmp -lgsl -O3
EXE=jmmOneDMC

henrun : henry2
	./${EXE} 100 0.5 0.5 -1 1000 lj 0.05 2 10 10 0.01 100000 100 10 5 0.01 100000 100 125 > out.out

henrun2 : henry2
	./${EXE} 1000 0.25 0.75 5 10000 lj 0.05 2 10 10 0.01 100000 100 10 5 0.01 100000 100 125 > out.out

henry2 : 
	${CChen} ${CFLAGShen} Main.c -o ./${EXE}

macbook :
	${CCmac} ${CFLAGSmac} Main.c -o ./${EXE}

clean :
	rm -rf *.o

