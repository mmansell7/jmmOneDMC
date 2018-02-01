/*
 ============================================================================
 Name        : Main.c
 Author      : Matt Mansell, NC State University, Dept. of Chemical &
                 Biomolecular Engineering
 Version     : 0.0.1
 Copyright   : Your copyright notice
 Description : Main function for One-Dimensional NPT Simulations
 ============================================================================
 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "jmmMCState.h"

typedef double* dVECTOR;
typedef dVECTOR* dMATRIX;
typedef unsigned long int* luVECTOR;
typedef luVECTOR* luMATRIX;
typedef int* iVECTOR;
typedef iVECTOR* iMATRIX;
typedef FILE** fVECTOR;

struct inputData {
	unsigned long int N,numSteps,cpi,tpi;
	double P,T,maxStep,maxdl;
	double * (*phi)(double);
};

struct mcState {
	dVECTOR r, rTrial,rij,rijTrial,e6ij,e6ijTrial,e12ij,e12ijTrial,vir6ij,vir6ijTrial,vir12ij,vir12ijTrial;
	double *lp,*Ep,*Virp;
};

int tid,bndchk;
unsigned long int stepNum,cpi,tpi,slcp,sltp;
unsigned long int N,numPairs,numSteps,sn; // sn = "step number"
unsigned long int ind,ii,jj,dj,nm,dAcc[2],vAcc[2];
double P,T,E,dE,ETrial,E6,dE6,E6Trial,E12,dE12,E12Trial,l,Vir,VirTrial,Vir6,Vir6Trial,Vir12,Vir12Trial,dVir,dVir6,dVir12;
double rij1,rij3,rij6,rij7,rij12,rij13;
double delta,maxStep,maxdl,rndm,md,dl,lRat1,lRat3,lRat6,lRat7,lRat12,lRat13,ETrial,VirTrial,ran,bf,lA,lSA,EA,ESA,VirA,VirSA;
dVECTOR r,rTrial,rij,rijTrial,e12ij,e12ijTrial,e6ij,e6ijTrial,ij,vir6ij,vir6ijTrial,vir12ij,vir12ijTrial;
luVECTOR iii,jjj;
luMATRIX indind;
FILE *cf, *tf;
char *verStr = "1.0.0";
char *verDateStr = "01/19/2018";
struct MCState *mcs;

// For g_i(x) and rho(x), I need...
int gns;
unsigned long int gpi,rhopi,gnb,rhonb,slg,slrho;
unsigned long int seed;
double rbw,gsw,gbw;
FILE *rhof;
fVECTOR gf;
iVECTOR rhol,rhoA;
dVECTOR rhoM;
iMATRIX gl,gA;
dMATRIX gM;
long int gs11,gs12;
// 0. Done. Read in rhobw,rhonb,rhopi,gsw,gns,gbw,gnb,gpi
// 1. Done. Allocate vectors and matrices, and set initial values
// 2. Done. For accepted displacements, accumulate
// 3. Done. For rejected displacements, accumulate
// 4. Done. For accepted volume changes, accumulate
// 5. Done. For rejected volume changes, accumulate
void fgrho();  // Full update of rho(x) accumulator and g(x) accumulator
// 6. For rhopi, print rho
void printRho();
// 7. For gpi, print g
void printG();

double * phiLJ(double *d, void *params);
double * (*phi)(double *d, void *params);
void fad(unsigned long int N, unsigned long int np, dVECTOR r, dVECTOR rTrial, dVECTOR rij, dVECTOR rijTrial, dVECTOR e6ij, dVECTOR e6ijTrial, dVECTOR e12ij, dVECTOR e12ijTrial,
		dVECTOR vir6ij, dVECTOR vir6ijTrial, dVECTOR vir12ij, dVECTOR vir12ijTrial, double *lp,double *Ep, double *Virp, unsigned long int nm,
		double rn);  // Full attempt, displacement
void qad();  // Quick attempt at displacement
void qav();  // Quick attempt at volume change
void fav(dVECTOR r, dVECTOR rTrial, dVECTOR rij, dVECTOR rijTrial, dVECTOR e6ij, dVECTOR e6ijTrial, dVECTOR e12ij, dVECTOR e12ijTrial,
		dVECTOR vir6ij, dVECTOR vir6ijTrial, dVECTOR vir12ij, dVECTOR vir12ijTrial, double *lp,double *Ep, double *Virp, double rn);  // Full attempt, volume change
void printCoords(FILE *cf, unsigned long int sn, dVECTOR r);
void printThermo(FILE *tf, unsigned long int *sltpp, unsigned long int sn, double *lAp, double *lSAp, double *EAp, double *ESAp, double *VirAp, double *VirSAp);
struct MCInput readInput(int argc, char *argv[],unsigned long int *N,double *P,double *T,unsigned long int *ns,double * (*phi)(double *,void *),double *ms,
		double *mdl,unsigned long int *cpi,unsigned long int *tpi);
void setup(FILE **cfp, FILE **tfp, double *Ep, double *Virp, double *lp, gsl_rng **rg, unsigned long int *slcp, unsigned long int *sltp,
		double *lA, double *lSA,	double *EA, double *ESA, double *VirA, double *VirSA, dVECTOR *r, dVECTOR *rTrial, dVECTOR *rij,
		dVECTOR *rijTrial, dVECTOR *e6ij, dVECTOR *e6ijTrial, dVECTOR *e12ij, dVECTOR *e12ijTrial, dVECTOR *vir12ij, dVECTOR *vir12ijTrial,
		dVECTOR *vir6ij, dVECTOR *vir6ijTrial,unsigned long int seed);
static gsl_rng * rangen;
void fagr(FILE *gf, FILE *rhof, struct mcState *mcsp);
void qagrho(int tid);  // Quick calculation of rho(x) accumulator and g(x) accumulator after displacement
void ugrho();  // Update accumulators


int main (int argc, char *argv[]) {

//	struct mcState mcs;
	printf("########################################\n");
	printf("#       Matt Mansell's MC Code         #\n");
	printf("#      Version %s : %s      #\n",verStr,verDateStr);
	printf("########################################\n\n");
	struct MCInput inp = readInput(argc, argv, &N,&P,&T,&numSteps,phi,&maxStep,&maxdl,&cpi,&tpi);
	printf("N: %lu\nP: %.5G\nT: %.5G\nnumSteps: %lu\nmaxStep: %.5G\nmax vol change: %.5G\n"N,P,T,numSteps,maxStep,maxdl);
	printf("Configuration print interval: %lu\nThermo print interval: %lu\nDensity bin width: %.5G\n",cpi,tpi,rbw);
	printf("Number of density bins: %lu\nDensity print interval: %lu\ng(x) (or two-particle density) segment width: %.5G\n",rnb,rhopi,gsw);
	printf("Number of g(x) segments: %lu\ng(x) bin width: %.5G\nNumber of g(x) bins: %lu\n",gns,gbw,gnb);
	printf("g(x) print interval: %lu\n\nSeed: %lu\n\n",gpip,seed);
	fflush(stdout);
	setup(&cf,&tf,&E,&Vir,&l,&rangen,&slcp,&sltp,&lA,&lSA,&EA,&ESA,&VirA,&VirSA,&r,&rTrial,&rij,&rijTrial,&e6ij,
			&e6ijTrial,&e12ij,&e12ijTrial,&vir12ij,&vir12ijTrial,&vir6ij,&vir6ijTrial,seed);
	mcs = setupMCS(inp);
	fprintf(stdout,"Setup completed\n");
	sn = 0;

#pragma omp parallel default (shared) \
private(tid,rij1,rij3,rij6,rij7,rij12,rij13) \
shared(mcs,E,dE6,dE12,ETrial,Vir,VirTrial,iii,jjj,l,lRat1,r,rij,e6ij,e12ij,vir6ij,vir12ij,lA,lSA,EA,ESA,nm,rndm)
{
	tid = omp_get_thread_num();

	#pragma omp single
	{
		printf("Running with %d threads...\n\n",omp_get_num_threads());
		mcsPrintStep(mcs);
		//printf("Step: %lu...\n",sn);
		fflush(stdout);
	}
	fad(N,numPairs,r,rTrial,rij,rijTrial,e6ij,e6ijTrial,e12ij,e12ijTrial,vir6ij,vir6ijTrial,vir12ij,vir12ijTrial,&l,&E,&Vir, 5, 0.5);
	unsigned long int ddd = 5;
	double nnn = 0.5;
	mcs_fad(mcs,&ddd,&nnn);
	
	#pragma omp single
	{
		printCoords(cf,sn,r);
	}
	#pragma omp single
	{
		slrho = -1;
		printRho();
	}
	
	slg = -1;
	printG();
	
	// ############################### MC LOOP #######################################
	while (sn < numSteps) {
		#pragma omp barrier
		#pragma omp single
		{
		sn++;
		if (sn % 1000 == 0) {
			printf("Step: %lu...\n",sn);
			//fflush(stdout);
		}
		nm = gsl_rng_uniform_int(rangen,N+1);
		rndm = gsl_rng_uniform(rangen);
		}

		if (nm < N) {
			qad();
			//fprintf(rhof,"nm: %lu, md: %.5G\n",nm,md);
			//fad(N,numPairs,r,rTrial,rij,rijTrial,e6ij,e6ijTrial,e12ij,e12ijTrial,vir6ij,vir6ijTrial,vir12ij,vir12ijTrial,&l,&E,&Vir,nm,rndm);
			}
		else {
			//fav(r,rTrial,rij,rijTrial,e6ij,e6ijTrial,e12ij,e12ijTrial,vir6ij,vir6ijTrial,vir12ij,vir12ijTrial,&l,&E,&Vir,   rndm);
			qav();
		}

		//fgrho();
		//printf("Entering sections\n");
		//fflush(stdout);
		#pragma omp barrier
		#pragma omp sections private(ii,jj,dj,ind) nowait
		{
			#pragma omp section
			{
				if (sn % cpi == 0) {
					//printf("Printing coordinates...\n");
					printCoords(cf,sn,r);
					//printf("Done printing coordinates...\n");
				}
			}

			#pragma omp section
			{
				lA = lA + l;
				lSA = lSA + l*l;
				EA = EA + E;
				ESA = ESA + E*E;
				VirA = VirA + Vir;
				VirSA = VirSA + Vir*Vir;
				//printf("\nE,EA = %lf,%lf\n",E,EA);
				//fflush(stdout);
				//ind = 0;
				//for (ii = 0; ii < N - 1; ii++) {
				//	for (jj = ii + 1; jj < N; jj++) {
						//fprintf(tf,"indind[%lu][%lu] = %lu\n",ii,jj,indind[ii][jj]);
						//ind++;
						//fprintf(tf,"indind[%lu][%lu] = %lu\n",ii,jj,indind[ii][jj]);
						//ind++;
				//	}
				//}
				if (sn % tpi == 0) {
					//printf("Printing thermodynamic data...\n");
					//for (ii = 0; ii < N - 1; ii++) {
					//	for (jj = ii + 1; jj < N; jj++) {
					//		dj = jj - ii - 1;
					//		ind = indind[ii][dj];
					//		ind = indind[ii][dj];
					//		fprintf(tf,"%lf ",e12ij[ind]);
						//,e6ij[ind],e12ij[ind]);
							//fflush(tf);
							//fprintf(tf,"tid: %d...ii,jj,dj,ind = %lu,%lu,%lu,%lu...rij[%lu] = %lf\n",tid,ii,jj,dj,ind,ind,rij[ind]);
							//fflush(tf);
					//	}
					//	fprintf(tf,"\n");
					//}
					//fprintf(tf,"\nE,EA = %lf,%lf\n",E,EA);
					//fflush(tf);
					printThermo(tf,&sltp,sn,&lA,&lSA,&EA,&ESA,&VirA,&VirSA);
					//printf("Done printing thermo data...\n");
				}
			}

			#pragma omp section
			{
				if (sn % rhopi == 0) {
					//printf("Printing rho...\n");
					printRho();
					//printf("Done printing rho...\n");
				}
			}
		}
			
		if (sn % gpi == 0) {
			//#pragma omp master
			//{
			//printf("Printing g(x)...\n");
			//}

			printG();
			//#pragma omp single
			//{
			//printf("Done printing g(x)\n");
			//}
		}
		
		if (sn % 10000 == 0) {
			unsigned long int iq;
			//ind = 0;
			ETrial = 0;
			#pragma omp for private(iq,ii,jj,rij1,rij3,rij6,rij12) reduction(+:ETrial)
			//for (ii = 0; ii < N - 1; ii++) {
				//for (jj = ii + 1; jj < N; jj++) {
				for (iq = 0; iq < numPairs; iq++) {
					//fprintf(stdout,"iq = %lu",iq);
					//fflush(stdout);
					ii = iii[iq];
					jj = jjj[iq];
					
					rijTrial[iq] = r[jj] - r[ii];
					rij1 = rijTrial[iq];
					rij3 = rij1*rij1*rij1;
					rij6 = 1/(rij3*rij3);
					rij12 = rij6*rij6;
					e12ijTrial[iq] = 4*rij12;
					e6ijTrial[iq] = 4*rij6;
					ETrial += e12ijTrial[iq] - e6ijTrial[iq];
					//fprintf(stdout,"iq = %lu\n",iq);
					//fflush(stdout);
					//iq++;
				//}
			}
			#pragma omp single
			{
			if (fabs(ETrial - E) > 0.001) {
				fprintf(stdout,"\n\n\nEnergy discrepancy!!! ######### !!!!!!!!!!! $$$$$$$$$$ !!!!!!!!\nE = %.8G\nETest = %.8G\nResetting energy to that of full calculation (E = %.8G)\n\n",E,ETrial,E);
				E = ETrial;
				/*ind = 0;
				for (ii = 0; ii < N; ii++) {
					fprintf(stdout,"%lu %.8G ",ii,r[ii]);
					for (jj = ii + 1; jj < N; jj++) {
						dj = jj - ii - 1;
						if (fabs(rijTrial[ind] - rij[ind]) > 0.001) {
							fprintf(stdout,"%lu,%lu : %.8G ",ii,jj,rijTrial[ind]-rij[ind]);
						}
						ind++;
					}
					fprintf(stdout,"\n");
				} */ 
				
			}  
			else {
				fprintf(stdout,"####### Energy was just verified ########  E = %.8G ######## ETest = %.8G #######\n",E,ETrial);
			}
			}
		}
		//printf("gA: %p %p gA[0]: %p %p  gA[%u]: %p %p\n",gA,&(gA[0]),gA[0],&(gA[0][0]),gns-1,gA[gns-1],&(gA[gns-1][0]));
		#pragma omp single
		{
		//printf("\n");
		}
		fflush(stdout);

	}     //  ############################### END MC LOOP #######################################

#pragma omp barrier

}

	printf("\nPROGRAM COMPLETED SUCCESSFULLY!\n");
	printf("\nE = %.8G\n",E);
	printf("dAcc dRej vAcc vRej\n%lu %lu %lu %lu\n",dAcc[0],dAcc[1],vAcc[0],vAcc[1]);
	return 0;

}


double * phiLJ(double *d, void *params) {
	double rij1,rij3,rij6,rij7,rij12,rij13,phi6,phi12,vir6,vir12;
	double *phi = (double *) malloc(4*sizeof(double));
	//printf("Check inside ####### d = %lf\n",d);

	rij1 = *d;
	rij3 = rij1*rij1*rij1;
	rij6 = 1/(rij3*rij3);
	rij7 = rij6/rij1;
	rij12 = rij6*rij6;
	rij13 = rij12/rij1;

	phi6 = 4*rij6;
	phi12 = 4*rij12;
	vir6 = 24/l*rij6;
	vir12 = 48/l*rij12;
	phi[0] = phi6;
	phi[1] = phi12;
	phi[2] = vir6;
	phi[3] = vir12;

	return phi;
}

void fad(unsigned long int N, unsigned long int np, dVECTOR r, dVECTOR rTrial, dVECTOR rij, dVECTOR rijTrial,
		dVECTOR e6ij, dVECTOR e6ijTrial, dVECTOR e12ij, dVECTOR e12ijTrial, dVECTOR vir6ij, dVECTOR vir6ijTrial, dVECTOR vir12ij,
		dVECTOR vir12ijTrial, double *lp,double *Ep, double *Virp, unsigned long int nm, double rn) {
	// Full attempt, displacement
	unsigned long int ii,jj,ind;
	double rij1,rij3,rij6,rij7,rij12,rij13;

	#pragma omp single
	{
	for (ind = 0; ind < N; ind++) {
		rTrial[ind] = r[ind];
	}
	for (ind = 0; ind < np; ind++) {
		rijTrial[ind] = rij[ind];
	}

	md = (rn - 0.5)*2*maxStep;
	//fprintf(stdout,"Attempting to move particle %lu by %lf...",nm,md);
	rTrial[nm] = r[nm] + md;
	E6Trial = 0;
	E12Trial = 0;
	ETrial = 0;
	Vir6Trial = 0;
	Vir12Trial = 0;
	VirTrial = 0;
//	printf("Particle to move: %lu,  Move distance: %.5f...",nm,md);
	if (fabs(rTrial[nm]) > (l/2.0)) {
		bndchk = 1;
		//printf("Wall collision. Displacement rejected\n");
		dAcc[1]++;
	}
	else {
		bndchk = 0;
	}
	fflush(stdout);
	}


	if (bndchk == 0) {
		#pragma omp for private(ii,jj,ind,rij1,rij3,rij6,rij7,rij12,rij13) reduction(+:E6Trial,E12Trial,Vir6Trial,Vir12Trial)
		for (ind = 0; ind < np; ind++) {
			ii            = iii[ind];
			jj            = jjj[ind];
			rijTrial[ind] = rTrial[jj] - rTrial[ii];
//			printf("ii,jj,rij = %lu,%lu, %.5G\n",ii,jj,rijTrial[ind]);
			rij1          = rijTrial[ind];
			rij3          = rij1*rij1*rij1;
			rij6          = 1/(rij3*rij3);
			rij7          = rij6/rij1;
			rij12         = rij6*rij6;
			rij13         = rij12/rij1;

			e6ijTrial[ind]     = 4*rij6;
			e12ijTrial[ind]    = 4*rij12;
			vir6ijTrial[ind]   = 24/l*rij6;
			vir12ijTrial[ind]  = 48/l*rij12;
			
			E6Trial     += e6ijTrial[ind];
			E12Trial    += e12ijTrial[ind];
			Vir6Trial   += vir6ijTrial[ind];
			Vir12Trial  += vir12ijTrial[ind];
//			printf("ind,ii,jj,rijTrial,e6Trial,e12Trial,ETrial = %lu,%lu,%lu, %.5G, %.5G, %.5G, %.5G\n",ind,ii,jj,rijTrial[ind],e6ijTrial[ind],e12ijTrial[ind],ETrial);

		}
//		printf("Done with for loop\n");

		#pragma omp single
		{
		ETrial = E12Trial - E6Trial;

		if (ETrial > (*Ep) ) {
//			printf("ETrial: %lf...*Ep: %.5G...",ETrial,*Ep);
			ran = gsl_rng_uniform(rangen);
			bf = exp((*Ep-ETrial)/T);
		}

		if (ETrial < *Ep || bf > ran) {
			//printf("Displacement accepted...E = %.8G\n",ETrial);
			dAcc[0]++;
			*Ep   = ETrial;
			E6    = E6Trial;
			E12   = E12Trial;
			Vir6  = Vir6Trial;
			Vir12 = Vir12Trial;
			*Virp = N*T/l + Vir12 - Vir6;
			r[nm] = rTrial[nm];
			fflush(stdout);
			for (ind = 0; ind < np; ind++) {
				rij[ind]     = rijTrial[ind];
				e6ij[ind]    = e6ijTrial[ind];
				e12ij[ind]   = e12ijTrial[ind];
				vir6ij[ind]  = vir6ijTrial[ind];
				vir12ij[ind] = vir12ijTrial[ind];
			}
		}
		else {
			//printf("Displacement rejected\n");
			dAcc[1]++;
		}
		//printf("Check 9...Thread no.: %d\n",omp_get_thread_num());
		}
	}
	//printf("Exiting fad...Thread no.: %d\n",omp_get_thread_num());
	fflush(stdout);
	#pragma omp barrier
}

void fav(dVECTOR r, dVECTOR rTrial, dVECTOR rij, dVECTOR rijTrial, dVECTOR e6ij, dVECTOR e6ijTrial, dVECTOR e12ij, dVECTOR e12ijTrial,
		dVECTOR vir6ij, dVECTOR vir6ijTrial, dVECTOR vir12ij, dVECTOR vir12ijTrial, double *lp,double *Ep, double *Virp, double rn) {

	unsigned long int ii,jj,ind;
	double rij1a,rij3a,rij6a,rij7a,rij12a,rij13a;

	#pragma omp single
	{
	bndchk = 0;
	ETrial = 0;
	VirTrial = 0;
	dl = (rn-0.5)*2*maxdl;
//	printf("Volume change attempt. Linitial = %.5G, dL = %.5G...",*lp,dl);
	lRat1 = (*lp+dl)/(*lp);
//	printf("dl = %.5G, lRat1 = %.5G\n",dl,lRat1);
	for (ind = 0; ind < N; ind++) {
		rTrial[ind] = r[ind]*lRat1;
	}

	}

	#pragma omp for private(ii,jj,ind,rij1a,rij3a,rij6a,rij7a,rij12a,rij13a) reduction (+:ETrial,VirTrial)
	for (ind = 0; ind < numPairs; ind++) {
		ii            = iii[ind];
		jj            = jjj[ind];

//		printf("ii,jj,rij,= %lu,%lu, %.5G\n",ii,jj,rij[ind]);

		rij1a          = rTrial[jj] - rTrial[ii];

//		printf("ii,jj,rij,= %lu,%lu, %.5G\n",ii,jj,rij[ind]);

		rij3a          = rij1a*rij1a*rij1a;
		rij6a          = 1/(rij3a*rij3a);
		rij7a          = rij6a/rij1a;
		rij12a         = rij6a*rij6a;
		rij13a         = rij12a/rij1a;

		e6ijTrial[ind]     = 4*rij6a;
		e12ijTrial[ind]    = 4*rij12a;
		vir6ijTrial[ind]     = 24/l*rij6a;
		vir12ijTrial[ind]    = 48/l*rij12a;
//		printf("ii,jj,rij1,rij6,e6 = %lu,%lu, %.5G, %.5G, %.5G\n",ii,jj,rij[ind],rij6a,e6ij[ind]);

		ETrial            += e12ijTrial[ind] - e6ijTrial[ind];
		VirTrial          += vir12ijTrial[ind] - vir6ijTrial[ind];
	}

	#pragma omp single
	{

//	printf("ETrial: %lf...",ETrial);
	bf = exp(-(ETrial - *Ep + P*dl)/T)*pow(lRat1,N);

	if (bf < 1.0) {
		ran = gsl_rng_uniform(rangen);
	}

	if (bf >= 1.0 || bf > ran) {
//		printf("Volume change accepted...");
		vAcc[0]++;
		*lp = *lp + dl;
		*Ep = ETrial;
		*Virp = VirTrial;
		for (ind = 0; ind < N; ind++) {
			r[ind] = rTrial[ind];
		}
		for (ind = 0; ind < numPairs; ind++) {
			rij[ind] = rijTrial[ind];
			e6ij[ind] = e6ijTrial[ind];
			e12ij[ind] = e12ijTrial[ind];
			vir6ij[ind] = vir6ijTrial[ind];
			vir12ij[ind] = vir12ijTrial[ind];
		}
//		printf("Volume now: %.5G\n",*lp);
	}
	else {
		vAcc[1]++;
//		printf("Volume change rejected\n");
	}
	fflush(stdout);
	}
	#pragma omp barrier
}

void printCoords(FILE *cf, unsigned long int sn, dVECTOR r) {
	fprintf(cf,"%lu\nStep no.: %lu  Box length: %.5f\n",N,sn,l);
	for(ind = 0; ind < N; ind++) {
		fprintf(cf,"%lu  0.0  0.0  %.8G\n",ind+1,r[ind]);
	}
	fflush(cf);
}

void printThermo(FILE *tf, unsigned long int *sltpp, unsigned long int sn, double *lAp, double *lSAp, double *EAp, double *ESAp, double *VirAp, double *VirSAp) {
	double lM,lSM,EM,ESM,VirM,VirSM;
	unsigned long int ss;

	ss = sn - *sltpp;

	lM       = *lAp/ss;
	lSM      = *lSAp/ss;
	EM       = *EAp/ss;
	ESM      = *ESAp/ss;
	VirM     = *VirAp/ss;
	VirSM    = *VirSAp/ss;
//	printf("Check 3...sn: %lu...EM: %.5G...ESM: %.5G...lM: %.5G...lSM: %.5G\n",sn,EM,ESM,lM,lSM);
	
	fprintf(tf,"%lu\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\n",sn,EM,ESM,lM,lSM,VirM,VirSM);
	fflush(tf);
	*lAp     = 0;
	*lSAp    = 0;
	*EAp     = 0;
	*ESAp    = 0;
	*VirAp   = 0;
	*VirSAp  = 0;
	*sltpp   = sn;

}

struct MCInput readInput(int argc, char *argv[],unsigned long int *N,double *P,double *T,unsigned long int *ns,double * (*phi)(double*,void*),double *ms,
		double *mdl,unsigned long int *cpi,unsigned long int *tpi) {

	/*	Read input arguments:
			1.  Number of Particles
			2.  Pressure
			3.  Temperature
			4.  Number of trial steps
			5.  Potential type
			6.  Maximum particle displacement
			7.  Maximum box length change
			8.  Configuration print interval
			9.  Thermo print interval
			10. Rho bin width
			11. Number of bins for rho
			12. Rho print interval
			13. g(x) segment width
			14. Number of segments for g(x)
			15. g(x) bin width
			16. Number of bins for g(x)
			17. g(x) print interval
			18. Seed (suggested use: date +'%s')

		Identified future extensions:

	*/
	
        struct MCInput inp;
       
	*N = (unsigned long int) strtol(argv[1],NULL,10);
	*P = strtod(argv[2],NULL);
	*T = strtod(argv[3],NULL);
	*ns = (unsigned long int) strtol(argv[4],NULL,10);
	if (strcmp(argv[5],"lj") == 0) {
		printf("Using Lennard-Jones potential\n");
		phi = phiLJ;
	}
	else {
		printf("Unknown potential type requested (%s). Using Lennard-Jones potential.\n",argv[5]);
		phi = phiLJ;
	}

	*ms      = (double) strtod(argv[6],NULL);
	*mdl     = (double) strtod(argv[7],NULL);
	*cpi     = (unsigned long int) strtol(argv[8],NULL,10);
	*tpi     = (unsigned long int) strtol(argv[9],NULL,10);

	rbw      = (double) strtod(argv[10],NULL);
	rhonb    = (unsigned long int) strtol(argv[11],NULL,10);
	rhopi    = (unsigned long int) strtol(argv[12],NULL,10);

	gsw      = (double) strtod(argv[13],NULL);
	gns      = (int) strtol(argv[14],NULL,10);
	gbw      = (double) strtod(argv[15],NULL);
	gnb      = (unsigned long int) strtol(argv[16],NULL,10);
	gpi      = (unsigned long int) strtol(argv[17],NULL,10);
	seed     = (unsigned long int) strtol(argv[18],NULL,10);
        
	inp.N = *N;
	inp.P = *P;
	inp.T = *T;
	inp.ns = *ns;
	inp.phi = phi;
	inp.gpi = gpi;
	inp.maxStep = *ms;
	inp.maxdl = *mdl;
	inp.cpi = *cpi;
	inp.tpi = *tpi;
	inp.rbw = rbw;
	inp.rhonb = rhonb;
	inp.rhopi = rhopi;
	inp.gsw = gsw;
	inp.gns = gns;
	inp.gbw = gbw;
	inp.gnb = gnb;
	
        return inp;
}

void setup(FILE **cfp, FILE **tfp, double *Ep, double *Virp, double *lp, gsl_rng **rg, unsigned long int *slcp, unsigned long int *sltp,
		double *lA, double *lSA,	double *EA, double *ESA, double *VirA, double *VirSA, dVECTOR *r, dVECTOR *rTrial, dVECTOR *rij,
		dVECTOR *rijTrial, dVECTOR *e6ij, dVECTOR *e6ijTrial, dVECTOR *e12ij, dVECTOR *e12ijTrial, dVECTOR *vir12ij, dVECTOR *vir12ijTrial,
		dVECTOR *vir6ij, dVECTOR *vir6ijTrial,unsigned long int seed) {
	//struct MCState *mcs;

	unsigned long int ii,jj,dj,ind;
	char gfstr[20];
	
	*rg = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(*rg,seed);
	
	*cfp = fopen("config.dat","w");
	*tfp = fopen("thermo.dat","w");
	fprintf(*tfp,"Step  Energy  Energy^2    l     l^2     Virial  Virial^2\n");

	*slcp = 0;
	*sltp = 0;
	*lA = 0;
	*lSA = 0;
	*EA = 0;
	*ESA = 0;
	*VirA = 0;
	*VirSA = 0;

	dAcc[0] = 0;
	dAcc[1] = 0;
	vAcc[0] = 0;
	vAcc[1] = 0;

	numPairs       = (unsigned long int) ((double) (N-1)/2*N);
	*Ep            = 10E10;
	*lp            = N;
	*Virp          = 0;

	*r             = (dVECTOR) malloc(N*sizeof(double));
	*rij           = (dVECTOR) malloc(numPairs*sizeof(double));
	*e12ij         = (dVECTOR) malloc(numPairs*sizeof(double));
	*e6ij          = (dVECTOR) malloc(numPairs*sizeof(double));
	*vir12ij         = (dVECTOR) malloc(numPairs*sizeof(double));
	*vir6ij          = (dVECTOR) malloc(numPairs*sizeof(double));

	*rTrial        = (dVECTOR) malloc(N*sizeof(double));
	*rijTrial      = (dVECTOR) malloc(numPairs*sizeof(double));
	*e12ijTrial    = (dVECTOR) malloc(numPairs*sizeof(double));
	*e6ijTrial     = (dVECTOR) malloc(numPairs*sizeof(double));
	*vir12ijTrial    = (dVECTOR) malloc(numPairs*sizeof(double));
	*vir6ijTrial     = (dVECTOR) malloc(numPairs*sizeof(double));

	iii            = (luVECTOR) malloc(numPairs*sizeof(double));
	jjj            = (luVECTOR) malloc(numPairs*sizeof(double));
	indind         = (luMATRIX) malloc((N-1)*sizeof(luVECTOR));
	
	ind = 0;
	for (ii = 0; ii < N-1; ii++) {
		indind[ii] = (luVECTOR) malloc((N-1-ii)*sizeof(unsigned long int));
		for (jj = ii + 1; jj < N; jj++) {
			dj = jj - ii - 1;
			indind[ii][dj] = ind;
			//fprintf(stdout,"ii,jj,dj,indind[%lu][%lu] = %lu,%lu,%lu,%lu\n",ii,jj,ii,jj,dj,indind[ii][dj]);
			ind++;
		}
	}
	
	slrho = 0;
	rhol = (iVECTOR) malloc(rhonb*sizeof(int));
	rhoA = (iVECTOR) malloc(rhonb*sizeof(int));
	rhoM = (dVECTOR) malloc(rhonb*sizeof(double));
	rhof = fopen("rho.dat","w");
	for (ii = 0; ii < N; ii++) {
		rhol[ii] = 0;
		rhoA[ii] = 0;
		rhoM[ii] = 0;
	}
	slg = 0;
	gl = (iMATRIX) malloc(gns*sizeof(iVECTOR));
	gA = (iMATRIX) malloc(gns*sizeof(iVECTOR));
	//printf("Just allocated gA at %p %p\n",gA,&gA[0]);
	gM = (dMATRIX) malloc(gns*sizeof(dVECTOR));
	gf = (fVECTOR) malloc(gns*sizeof(FILE *));
	for (ii = 0; ii < gns; ii++) {
		gl[ii] = (iVECTOR) malloc(gnb*sizeof(int));
		gA[ii] = (iVECTOR) malloc(gnb*sizeof(int));
		//printf("Just allocated gA[%lu] at %p %p\n",ii,gA[ii],(&gA[ii][0]));
		gM[ii] = (dVECTOR) malloc(gnb*sizeof(double));
		sprintf(gfstr,"g%lu.dat",ii);
		gf[ii] = fopen(gfstr,"w");
		for (jj = 0; jj < gnb; jj++) {
			gl[ii][jj] = 0;
			gA[ii][jj] = 0;
			gM[ii][jj] = 0;
		}
	}
	ind = 0;
	for (ii = 0; ii < N; ii++) {
		(*r)[ii]        = -l/2 + (ii+0.5)*(*lp/N);
		for (jj = ii+1; jj < N; jj++) {
			iii[ind] = ii;
			jjj[ind] = jj;
			ind++;
		}
	}
	for (ind = 0; ind < numPairs; ind++) {
		(*rij)[ind] = (*r)[jjj[ind]] - (*r)[iii[ind]];
	}

	fgrho();
	ugrho();
	fflush(stdout);
	
	//return mcs;
}

void fgrho() {

	unsigned long int ii,jj,ind;
	long int rb,gs1,gs2,gb;
	
	#pragma omp for
	for (rb = 0; rb < rhonb; rb++) {
		rhol[rb] = 0;
	}

	#pragma omp for
	for(ii = 0; ii < N; ii++) {
		//printf("ii = %lu\n",ii);
		rb = (long int) floor(r[ii]/rbw + rhonb/2.0);
		if (rb >= 0 && rb < rhonb) {
			#pragma omp atomic
			rhol[rb]++;
		}
	}
	
	//printf("Inside fgrho...Thread no.: %d...gns,gnb = %d, %lu\n",omp_get_thread_num(),gns,gnb);
	fflush(stdout);
	#pragma omp barrier
	// #pragma omp single
	// {
	
	#pragma omp for
	for (gs1 = 0; gs1 < gns; gs1++) {
		for (gb = 0; gb < gnb; gb++) {
			gl[gs1][gb] = 0;
		}
	}
	
	#pragma omp for private(ind,ii,jj,gs1,gs2,gb)
	//printf("gA: %p %p gA[0]: %p %p  gA[%u]: %p %p\n",gA,&(gA[0]),gA[0],&(gA[0][0]),gns-1,gA[gns-1],&(gA[gns-1][0]));
	for (ind = 0; ind < numPairs; ind++) {
		ii = iii[ind];
		jj = jjj[ind];
		gs1 = (long int) floor(r[ii]/gsw + gns/2.0);
		gs2 = (long int) floor(r[jj]/gsw + gns/2.0);
		gb = (long int) floor(rij[ind]/gbw);
		//printf("ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
		//fflush(stdout);
		if(gb < gnb) {
			if (gs1 >= 0 && gs1 < (long int) gns) {
				//printf("A ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
				//fflush(stdout);
				#pragma omp atomic
				gl[gs1][gb]++;
			}
			else {
				//printf("B ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
				//fflush(stdout);
				
			}
			
			if (gs2 >= 0 && gs2 < (long int) gns) {
				//printf("C ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
				//fflush(stdout);
				//printf("Accumulating to gA[%ld][%lu] at ",gs2,gb);
				//fflush(stdout);
				//printf("gA: %p  gA[%ld][%lu]: %p...",gA,gs2,gb,&(gA[gs2][gb]));
				//fflush(stdout);
				#pragma omp atomic
				gl[gs2][gb]++;
				//printf(" = %d\n",gA[gs2][gb]);
				//fflush(stdout);
			}
			else {
				//printf("D ");
				//printf("ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
				//fflush(stdout);
			}
		
		}
		else {
			//printf("E ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
			//fflush(stdout);
		}
		//printf("Next ind = %lu\n",ind+1);
	}
	//printf("Out of for loop\n");
}

void printRho() {
	unsigned long int rb,ns;
	//if (fprintf(rhof,"") < 0 ) {
	//	perror("Error printing to rhof.dat!\n");
	//}
	ns = (sn - slrho);
	//printf("printRho.  sn = %lu  rhof = %p  ns = %lu\n",sn,rhof,ns);
	//fflush(stdout);

	//printf("rhof new line\n");
	fprintf(rhof,"%lu",sn);
	//fflush(rhof);
	for (rb = 0; rb < rhonb; rb++) {
		//printf("## rhoM[%lu]: %p rhoA[%lu]: %p gA: %p %p gA[0]: %p %p  gA[%u]: %p %p\n",rb,&rhoM[rb],rb,&rhoA[rb],gA,&(gA[0]),gA[0],&(gA[0][0]),gns-1,gA[gns-1],&(gA[gns-1][0]));
		//fflush(stdout);
		//printf("rb = %lu\n",rb);
		//fflush(stdout);
		rhoM[rb] = (double) rhoA[rb]/ns/rbw;
		fprintf(rhof," %.8G",rhoM[rb]);
		//fflush(stdout);
		rhoA[rb] = 0;
	}
	fprintf(rhof,"\n");
	fflush(rhof);
	slrho = sn;
}

void printG() {
	int gs;
	unsigned long int gb,ns;
	
	#pragma omp for private(gs,gb,ns)
	for (gs = 0; gs < gns; gs++) {
		fprintf(gf[gs],"%lu",sn);
		ns = (sn - slg);

		for (gb = 0; gb < gnb; gb++) {
			gM[gs][gb] = (double) gA[gs][gb]/ns/gsw/gbw;
			fprintf(gf[gs]," %.8G",gM[gs][gb]);
			gA[gs][gb] = 0;
		}
		fprintf(gf[gs],"\n");
		fflush(gf[gs]);
	}
	slg = sn;
}


void qagrho(int tid) {
	
	unsigned long int ii,jj,rb,gs,gb,gb1,gb2;
	long int gs2;
	int aa,bb;
	//printf("Check 1...Thread %d...\n",tid);
	//fflush(stdout);
		
	#pragma omp single
	{
	long int rbn1,rbn2;
	rbn1 = (long int) floor((r[nm] - md)/rbw + rhonb/2.0);
	rbn2 = (long int) floor(r[nm]/rbw + rhonb/2.0);
	if (rbn1 >= 0 && rbn1 < rhonb) {
		rhol[rbn1]--;
	}
	if (rbn2 >= 0 && rbn2 < rhonb) {
		rhol[rbn2]++;
	}
	/*fprintf(rhof,"\n%lu: %.5G -> %.5G...%lu -> %lu\nrhol: ",nm,r[nm]-md,r[nm],rbn1,rbn2);
	for (rb = 0; rb < rhonb; rb++) {
		fprintf(rhof,"%d ",rhol[rb]);
	}
	fprintf(rhof,"\n");*/
	
	}
	
	#pragma omp single
	{	
	gs11 = (long int) floor((r[nm] - md)/gsw + gns/2.0);
	gs12 = (long int) floor(       r[nm]/gsw + gns/2.0);
	}
	
	aa = 0;
	bb = 0;
	if (gs11 >= 0 && gs11 < gns) aa += 1;
	if (gs12 >= 0 && gs12 < gns) aa += 2;

	#pragma omp for private(ii,jj,gs2,gb1,gb2) nowait
	for (ii = 0; ii < nm; ii++) {
//		printf("Check 2...Thread %d...ii %lu\n",tid,ii);
//		fflush(stdout);	
        	bb = aa;
	        gs2 = (long int) floor(r[ii]/gsw + gns/2.0);
		if (gs2 >= 0 && gs2 < gns) bb += 4;
		gb1 = (unsigned long int) floor(fabs((r[nm] - md - r[ii]))/gbw);
		if (gb1 >= 0 && gb1 < gnb) bb += 8;
		gb2 = (unsigned long int) floor(fabs((r[nm]      - r[ii]))/gbw);
		if (gb2 >= 0 && gb2 < gnb) bb += 16;

		//printf("Check 3...Thread %d...ii %lu...aa: %d...bb: %d\n",tid,ii,aa,bb);
		//fflush(stdout);	
		if (bb % 2 >= 1 && bb % 16 >= 8) {
			#pragma omp atomic
			gl[gs11][gb1]--;
		}
		if (bb % 4 >= 2 && bb % 32 >= 16) {
			#pragma omp atomic
			gl[gs12][gb2]++;
		}
		if (bb % 8 >= 4 && bb % 16 >= 8) {
			#pragma omp atomic
			gl[gs2][gb1]--;
		}
		if (bb % 8 >= 4 && bb % 32 >= 16) {
			#pragma omp atomic
			gl[gs2][gb2]++;
		}
//		printf("Check 4...Thread %d...ii %lu\n",tid,ii);
//		fflush(stdout);	
	}
	
	#pragma omp for private(ii,jj,gs2,gb1,gb2)
	for (ii = nm + 1; ii < N; ii++) {
        	bb = aa;
	        gs2 = (long int) floor(r[ii]/gsw + gns/2.0);
		if (gs2 >= 0 && gs2 < gns) bb += 4;
		gb1 = (unsigned long int) floor(fabs((r[nm] - md - r[ii]))/gbw);
		if (gb1 >= 0 && gb1 < gnb) bb += 8;
		gb2 = (unsigned long int) floor(fabs((r[nm]      - r[ii]))/gbw);
		if (gb2 >= 0 && gb2 < gnb) bb += 16;

		if (bb % 2 >= 1 && bb % 16 >= 8) {
			#pragma omp atomic
			gl[gs11][gb1]--;
		}
		if (bb % 4 >= 2 && bb % 32 >= 16) {
			#pragma omp atomic
			gl[gs12][gb2]++;
		}
		if (bb % 8 >= 4 && bb % 16 >= 8) {
			#pragma omp atomic
			gl[gs2][gb1]--;
		}
		if (bb % 8 >= 4 && bb % 32 >= 16) {
			#pragma omp atomic
			gl[gs2][gb2]++;
		}
	}
	
}


void qad() {

	unsigned long int ii,jj,dj,ind;
	//double d;
	
	md = (rndm - 0.5)*2*maxStep;
	
	#pragma omp single
	{
	rTrial[nm] = r[nm] + md;
	//fprintf(stdout,"Attempting to move particle %lu by %lf...",nm,md);
	}
	
	if (fabs(rTrial[nm]) > l/2.0) {
		#pragma omp single
		{
		dAcc[1]++;
		//printf("Wall collision. Move rejected...dRej = %lu",dAcc[1]);
		}
	}
	
	else {
		#pragma omp single
		{
		dE6 = 0;
		dE12 = 0;
		dVir6 = 0;
		dVir12 = 0;
		}
		
		//printf("nm = %lu...\n",nm);
		//fflush(stdout);	
		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) reduction(+:dE6,dE12,dVir6,dVir12) nowait
		for (ii = 0; ii < nm; ii++) {
		//	printf("ii = %lu\n",ii);
		//	fflush(stdout);	
			dj = nm - ii - 1;
			ind = indind[ii][dj];
			rijTrial[ind] = rij[ind] + md;
			rij1 = rijTrial[ind];
			rij3 = rij1*rij1*rij1;
			rij6 = 1/(rij3*rij3);
			rij7 = rij6/rij1;
			rij12 = rij6*rij6;
			rij13 = rij12/rij1;
			e6ijTrial[ind] = 4*rij6;
			e12ijTrial[ind] = 4*rij12;
			vir6ijTrial[ind] = 24/l*rij6;
			vir12ijTrial[ind] = 48/l*rij12;
			//printf("nm,ii,e6,e6Trial,e12,e12Trial = %lu,%lu, %lf, %lf, %lf, %lf\n",nm,ii,e6ij[ind],e6ijTrial[ind],e12ij[ind],e12ijTrial[ind]);
			dE6  = dE6  - e6ij[ind]  + e6ijTrial[ind];
			dE12 = dE12 - e12ij[ind] + e12ijTrial[ind];
			dVir6 = dVir6 - vir6ij[ind] + vir6ijTrial[ind];
			dVir12 = dVir12 - vir12ij[ind] + vir12ijTrial[ind];
		}
		
		//#pragma omp single
		//{
		//	if (nm == 19) {
		//		fprintf(stdout,"(Intermediate) dE6, dE12 = %.8G, %.8G\n",dE6,dE12);
		//	}
		//}
		
		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) reduction(+:dE6,dE12,dVir6,dVir12)
		for (ii = nm + 1; ii < N; ii++) {
		//	printf("ii = %lu\n",ii);
		//	fflush(stdout);	
			dj = ii - nm - 1;
			ind = indind[nm][dj];
			rijTrial[ind] = rij[ind] - md;
			rij1 = rijTrial[ind];
			rij3 = rij1*rij1*rij1;
			rij6 = 1/(rij3*rij3);
			rij7 = rij6/rij1;
			rij12 = rij6*rij6;
			rij13 = rij12/rij1;
			e6ijTrial[ind] = 4*rij6;
			e12ijTrial[ind] = 4*rij12;
			vir6ijTrial[ind] = 24/l*rij6;
			vir12ijTrial[ind] = 48/l*rij12;
			//printf("nm,ii,e6,e6Trial,e12,e12Trial = %lu,%lu, %lf, %lf, %lf, %lf\n",nm,ii,e6ij[ind],e6ijTrial[ind],e12ij[ind],e12ijTrial[ind]);
			dE6  = dE6  - e6ij[ind]  + e6ijTrial[ind];
			dE12 = dE12 - e12ij[ind] + e12ijTrial[ind];
			dVir6 = dVir6 - vir6ij[ind] + vir6ijTrial[ind];
			dVir12 = dVir12 - vir12ij[ind] + vir12ijTrial[ind];
		}
		
		//#pragma omp single
		//{
		//	if (nm == 19) {
		//		fprintf(stdout,"(End) dE6, dE12 = %.8G, %.8G\n",dE6,dE12);
		//	}
		//}
		#pragma omp barrier  // THIS BARRIER IS CRITICAL TO CALCULATE dE ACCURATELY (I DO NOT UNDERSTAND WHY THAT SHOULD BE THE CASE)
		
		#pragma omp single
		{
		dE = dE12 - dE6;
		//printf("dE = %.8G...",dE);
		//fflush(stdout);	
		}
		
		if (dE > 0 ) {
		#pragma omp single
		{
			//printf("ETrial: %lf...*Ep: %.5G...",ETrial,*Ep);
			ran = gsl_rng_uniform(rangen);
			bf = exp(-dE/T);
		}
		}
		
		if (dE <= 0 || bf > ran) {
			#pragma omp single
			{
			dAcc[0]++;
                        //printf("Boltzmann factor: %lf...random: %lf...",bf,ran);
			//printf("Displacement accepted...dAcc = %lu...",dAcc[0]);
			E      += dE;
			E6     += dE6;
			E12    += dE12;
			Vir6   += dVir6;
			Vir12  += dVir12;
			Vir    += dVir12 - dVir6;
			r[nm] = rTrial[nm];
			}
			
			#pragma omp for private(jj,dj,ind) nowait
			for (jj = 0; jj < nm; jj++) {
				//printf("jj = %lu\n",jj);
				//fflush(stdout);	
				dj = nm - jj - 1;
				ind = indind[jj][dj];
				//printf("1. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				rij[ind] = rijTrial[ind];
				//printf("2. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				e6ij[ind] = e6ijTrial[ind];
				//printf("3. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				e12ij[ind] = e12ijTrial[ind];
				//printf("4. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				vir6ij[ind] = vir6ijTrial[ind];
				//printf("5. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				vir12ij[ind] = vir12ijTrial[ind];
				//printf("6. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				//fprintf(tf,"jj,nm,dj,rij = %lu,%lu,%lu,%lf\n",jj,nm,dj,rij[ind]);
			}
			#pragma omp for private(jj,dj,ind)
			for (jj = nm + 1; jj < N; jj++) {
				//printf("7. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				dj = jj - nm - 1;
				//printf("8. nm,jj,dj,ind = %lu,%lu,%lu,%lu\n",nm,jj,dj,ind);
				//fflush(stdout);	
				ind = indind[nm][dj];
				//printf("9. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				rij[ind] = rijTrial[ind];
				//printf("10. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				e6ij[ind] = e6ijTrial[ind];
				//printf("11. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				e12ij[ind] = e12ijTrial[ind];
				//printf("12. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				vir6ij[ind] = vir6ijTrial[ind];
				//printf("13. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				vir12ij[ind] = vir12ijTrial[ind];
				//printf("14. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				//fprintf(tf,"jj,nm,dj,rij = %lu,%lu,%lu,%lf\n",jj,nm,dj,rij[ind]);
				fflush(tf);
			}
			
			qagrho(tid);
			
		/*#pragma omp single
		{
			//printf("Done accepting displacement\n");
			fflush(stdout);	
		}*/
		}
		else {
		#pragma omp single
		{
			dAcc[1]++;
                        //printf("Boltzmann factor: %lf...random: %lf...",bf,ran);
			//printf("Displacement rejected...dRej = %lu...",dAcc[1]);
		}
		}
		//printf("Check 9...Thread no.: %d\n",omp_get_thread_num());
		
		/*#pragma omp single
		{
			printf("E = %.8G\n",E);
		} */
	}
	
	ugrho();
	//printf("Exiting fad...Thread no.: %d\n",omp_get_thread_num());
	fflush(stdout);
	#pragma omp barrier
}




void qav() {
	
	unsigned long int ii,ind;
	double dQ;
	
	#pragma omp single
	{
	dl = (rndm-0.5)*2*maxdl;
	//fprintf(stdout,"Attempting to change volume by %lf...",dl);
	//printf("E6 = %.5G...E12 = %.5G...E = %.5G...",E6,E12,E);
	lRat1 = (l+dl)/l;
	lRat3 = lRat1*lRat1*lRat1;
	lRat6 = 1/(lRat3*lRat3);
	lRat12 = lRat6*lRat6;
	E6Trial = lRat6*E6;
	E12Trial = lRat12*E12;
	dE = E12Trial - E6Trial - E;
	//bf = exp(-(dE + P*dl)/T)*pow(lRat1,N);
	bf   = exp(-(dE + P*dl)/T + N*log(lRat1));
	//printf("E6Trial = %.5G...E12Trial = %.5G...dE = %.5G...Boltzmann factor: %lf",E6Trial,E12Trial,dE,bf);
	//printf("dE = %.5G...Boltzmann factor: %lf...",dE,bf);

	if (bf < 1.0 ) {
//		printf("ETrial: %lf...*Ep: %.5G...",ETrial,*Ep);
		ran = gsl_rng_uniform(rangen);
		//printf("dE: %lf...Boltzmann factor: %lf...",dE,bf);
	}
	}
	
	if (bf >= 1.0 || bf > ran) {
		#pragma omp single
		{
		//printf("random: %lf...",ran);
		vAcc[0]++;
		//printf("Volume change accepted...vAcc = %lu",vAcc[0]);
		E6  = E6Trial;
		E12 = E12Trial;
		E      = E12 - E6;
		
		l = l + dl;
		
		lRat7 = lRat6/lRat1;
		lRat13 = lRat12/lRat1;
		Vir6 = lRat7*Vir6;
		Vir12 = lRat13*Vir12;
		Vir    = N*T/l + Vir12 - Vir6;
		}
		
		#pragma omp single private(ii) nowait
		for (ii = 0; ii < N; ii++) {
			r[ii] = lRat1*r[ii];
		}
		
		#pragma omp for private(ind)
		for (ind = 0; ind < numPairs; ind++) { 
			rij[ind]     = lRat1*rij[ind];
			e6ij[ind]    = lRat6*e6ij[ind];
			e12ij[ind]   = lRat12*e12ij[ind];
			vir6ij[ind]  = lRat7*vir6ij[ind];
			vir12ij[ind] = lRat13*vir12ij[ind];
		}
		
/*		#pragma omp single
		{
		for (ii = 0; ii < N - 1; ii++) {
			for (jj = ii + 1; jj < N; jj++) {
				dj = jj - ii - 1;
				ind = indind[ii][dj];
				fprintf(stdout,"rij[%lu][%lu] = %lf\n",ii,jj,rij[ind]);
			}
		}
		} */
		fgrho();
	}
	
	else {
		#pragma omp single
		{
		//printf("random: %lf...",ran);
		vAcc[1]++;
		//printf("Volume change rejected...vRej = %lu",vAcc[1]);
		}
	}
	ugrho();
	fflush(stdout);
}



void ugrho() {
	
	unsigned long int rb,gb;
	unsigned int gs;
	
	#pragma omp for
	for (rb = 0; rb < rhonb; rb++) {
		//fprintf(rhof,"rhol[%lu] = %d ",rb,rhol[rb]);
		rhoA[rb] += rhol[rb];
	}
	/*#pragma omp single
	{
	fprintf(rhof,"rhoA: ");
	for (rb = 0; rb < rhonb; rb++) {
		fprintf(rhof,"%d ",rhoA[rb]);
	}
	fprintf(rhof,"\n");
	}*/
     
	//if (rbn1 != rbn2) {
		//fprintf(rhof,"md: %.5G, r[%lu]: %.5G -> %.5G,  rbn1 = %lu  rbn2 = %lu\n",md,nm,r[nm]-md,r[nm],rbn1,rbn2);
		//fprintf(rhof,"nm: %lu  rbn1 = %lu  rbn2 = %lu  rhol[rbn1] = %d   rhol[rbn2] = %d  rhoA[rbn1] = %d  rhoA[rbn2] = %d\n",nm,rbn1,rbn2,rhol[rbn1],rhol[rbn2],rhoA[rbn1],rhoA[rbn2]);
	//}

	#pragma omp for private(gs,gb)
	for (gs = 0; gs < gns; gs++) {
		for (gb = 0; gb < gnb; gb++) {
			gA[gs][gb] += gl[gs][gb];
		}
	}
}






