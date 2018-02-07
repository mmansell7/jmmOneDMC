#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "jmmMCState.h"

static double mcsETest,mcsETrial,mcsE6Trial,mcsE12Trial,mcsVirTrial,mcsVir6Trial,mcsVir12Trial;
static double mcsdE6,mcsdE12,mcsdE,mcsdVir6,mcsdVir12,mcsdVir,ran;
static double *rijTest,*e6ijTest,*e12ijTest;
static long int gs11,gs12;

struct MCState {
  
  // Use of typedefs to hide pointers is discouraged. Cease this practice.
  /*
  typedef double* dVECTOR;
  typedef dVECTOR* dMATRIX;
  typedef unsigned long int* luVECTOR;
  typedef luVECTOR* luMATRIX;
  typedef int* iVECTOR;
  typedef iVECTOR* iMATRIX;
  typedef FILE** fVECTOR;
  */


  unsigned long int N;             // Number of particles
  int nbn;                         // Number of neighbors with which each particle can
                                   //   interact.  -1 for no limit.
  unsigned long int sn;            // Step number
  unsigned long int numSteps;      // Number of steps (total steps requested)
  unsigned long int cpi;           // Configuration print interval
  unsigned long int tpi;           // Thermo print interval
  unsigned long int eci;           // Energy check interval
  unsigned long int mdai;          // Maximum displacement adjust interval
  unsigned long int mvai;          // Maximum volume change adjust interval
  unsigned long int slcp;          // Steps since last configuration print
  unsigned long int sltp;          // Steps since last thermo print
  unsigned long int numPairs;      // Number of pairs
  unsigned long int nm;            // Number of the particle selected for current trial
                                   //   nm = N means volume change
  unsigned long int gpi;           // g(r) print interval
  unsigned long int rhopi;         // Rho print interval
  unsigned long int gnb;           // Number of bins for g(r)
  unsigned long int rhonb;         // Number of bins for rho(r)
  unsigned long int slg;           // Steps since last g(r) print
  unsigned long int slrho;         // Steps since last rho(r) print
  unsigned long int seed;          // Seed value used whenever a random number generator
                                   //   is seeded or re-seeded
  unsigned long int dAcc[2];       // Array counting the trial moves accepted (dAcc[0]) and
                                   //   rejected (dAcc[1] )
  unsigned long int vAcc[2];       // Array counting the trial volume changes accepted
                                   //   (vAcc[0]) and rejected (vAcc[1])
  double P;                        // Pressure (reduced units)
  double T;                        // Temperature (reduced units)
  double maxStep;                  // Maximum trial displacement length
  double maxdl;                    // Maximum volume change trial
  double * (*phi)(double *rij,void *params);         // Pointer to the potential function
  double E;                        // Energy
  double dE;                       // Change in E for current trial
  double ETrial;                   // E for current trial
  double E6;                       // Total contribution (from all pairs) to energy from
                                   //   the ^(1/6) term in the Lennard-Jones potential
  double dE6;                      // Change in E6 for the current trial
  double E6Trial;                  // E6 for current trial
  double E12;                      // Total contribution (from all pairs) to energy from
                                   //   the ^(1/12) term in the Lennard-Jones potential
  double dE12;                     // Change in E12 for the current trial
  double E12Trial;                 // E12 for the current trial
  double l;                        // Simulation box length
  double Vir;                      // Virial
  double VirTrial;                 // Vir for current trial
  double Vir6;                     // Total contribution (from all pairs) to the virial
                                   //   from the ^(1/6) term in the Lennard-Jones potential
  double Vir6Trial;                // Vir6 for the current trial
  double Vir12;                    // Total contribution (from all pairs) to the virial
                                   //   from the ^(1/12) term in the Lennard-Jones potential
  double Vir12Trial;               // Vir12 for the current trial
  double dVir;                     // Change in Vir for the current trial
  double dVir6;                    // Change in Vir6 for the current trial
  double dVir12;                   // Change in Vir12 for the current trial

  double md;                       // Move distance for current trial displacement
  double rn;                       // Random number used to set md or dl
  double dl;                       // Length change for current trial displacement
  double bf;                       // Boltzmann factor for current trial
  double lA;                       // l accumulator
  double lSA;                      // l^2 accumulator
  double EA;                       // E accumulator
  double ESA;                      // E^2 accumulator
  double VirA;                     // Vir accumulator
  double VirSA;                    // Vir^2 accumulator
  double rbw;                      // rho(r) bin width
  double gsw;                      // g(r) segment width
  double gbw;                      // g(r) bin width
  
  int gns;                         // Number of segments for g(r)
  int bndchk;                      // Boundary check (0: no wall collision, 1: wall collision)
  
  double *r;                       // Vector of particle positions. Particles retain their
                                   //   numbers throughout the simulation. So, particle 0 is
                                   //   located at r[0], particle 1 at r[1], and so on up to
                                   //   particle N-1, located at r[N-1].
                                   
  double *rTrial;                  // r for current trial
                                   
  double *rij;                     // Vector of current pair separation distances pairs are 
                                   //   indexed sequentially by the number of the first particle,
                                   //   then the number of the second.  So, particle 0 and 
                                   //   particle 1 (pair 0-1) are indexed as 0 and are 
                                   //   separated by rij[0].  Pair 0-10 is separated by rij[9].
                                   //   And pair i-j is separated by rij[ i*(N-1) - (i*(i+1))/2
                                   //   + j - 1.
                                    
  double *rijTrial;                // rij for current trial
  double *e12ij;                   // Vector of pair contributions to E12
  double *e12ijTrial;              // e12ij for current trial
  double *e6ij;                    // Vector of pair contributions to E6
  double *e6ijTrial;               // e6ij for current trial 
  double *vir6ij;                  // Vector of pair contributions to Vir6
  double *vir6ijTrial;             // vir6ij for current trial
  double *vir12ij;                 // Vector of pair contributions to Vir12
  double *vir12ijTrial;            // vir12ij for current trial

  unsigned long int *iii;          // Vector of iii values for each pair.  iii[ind] is the index
                                   //   of the first particle in the pair.
  unsigned long int *jjj;          // As with iii, but gives the index of the second particle
  unsigned long int **indind;      // Inverse of iii and jjj -> indind[i][j] gives the index for
                                   //   the pair i-j

  FILE *cf;                        // Pointer to the configuration output file
  FILE *tf;                        // Pointer to the thermo output file
  FILE *rhof;                      // Pointer to the rho(r) output file
  FILE **gf;                       // Array of pointers to the g(r) files
  int *rhol;                       // Array of counts in each rho(r) bin after the last accepted
                                   //   trial.
  int *rhoA;                       // Array of accumulated counts in each rho bin since the start
                                   //   of the current rho(r) block
                                   
  double *rhoM;                    // Array describing the mean value (for the block) of rho(r) in
                                   //   each rho(r) bin
  int **gl;                        // Two-dimensional array with counts in each segment and bin of
                                   //   g(r) after the last accepted trial.
  int **gA;                        // Two-dimensional array with accumulated counts in each
                                   //   segment and bin of g(r) for the current g(r) block.
  double **gM;                     // Two-dimensional array with mean g(r) in each segment and bin
                                   //   for the current g(r) block
  gsl_rng * rangen;                // Random number generator
};

int printMCP(struct MCState *mcs1) {
	printf("Printing Monte Carlo parameters...\n");
        printf("N: %lu\nP: %.5G\nT: %.5G\n",mcs1->N,mcs1->P,mcs1->T);
	if ( mcs1->nbn > 0 ) {
		printf("Number of neighbors with which each particle can interact: %d\n",mcs1->nbn);
	}
	else {
		printf("No neighbor number limit.\n");
	}
	printf("numSteps: %lu\nmaxStep: %.5G\nmax vol change: %.5G\n",mcs1->numSteps,mcs1->maxStep,mcs1->maxdl);
        printf("Configuration print interval: %lu\nThermo print interval: %lu\nDensity bin width: %.5G\n",mcs1->cpi,mcs1->tpi,mcs1->rbw);
        printf("Number of density bins: %lu\nDensity print interval: %lu\ng(x) \
                  (or two-particle density) segment width: %.5G\n",mcs1->rhonb,
                  mcs1->rhopi,mcs1->gsw);
        printf("Number of g(x) segments: %d\ng(x) bin width: %.5G\nNumber of g(x) bins: %lu\n",mcs1->gns,mcs1->gbw,mcs1->gnb);
        printf("g(x) print interval: %lu\n\nSeed: %lu\n\n",mcs1->gpi,mcs1->seed);
        fflush(stdout);
	
	return 0;
}	

struct MCState * setupMCS(struct MCInput inp) {
        struct MCState *mcs = (struct MCState *) malloc(sizeof(struct MCState));

        unsigned long int ii,jj,dj,ind;
        char gfstr[20];
        
        mcs->N          = inp.N;
	mcs->P          = inp.P;
	mcs->T          = inp.T;
	mcs->nbn        = inp.nbn;
	mcs->sn         = 0;
        mcs->numSteps   = inp.ns;
        mcs->cpi        = inp.cpi;
        mcs->tpi        = inp.tpi;
        mcs->eci        = inp.eci;
        mcs->mdai       = inp.mdai;
        mcs->mvai       = inp.mvai;
	mcs->slcp       = -1;
	mcs->sltp       = -1;
	mcs->numPairs   = (unsigned long int) ((double) (mcs->N-1)/2*mcs->N);
	mcs->nm         = 0;
	mcs->gpi        = inp.gpi;
	mcs->rhopi      = inp.rhopi;
	mcs->gnb        = inp.gnb;
	mcs->rhonb      = inp.rhonb;
	mcs->maxStep    = inp.maxStep;
	mcs->maxdl      = inp.maxdl;
	mcs->phi        = inp.phi;
	mcs->rbw        = inp.rbw;
	mcs->gsw        = inp.gsw;
	mcs->gbw        = inp.gbw;
	mcs->gns        = inp.gns;
	mcs->seed       = inp.seed;

	printMCP(mcs);
	
	mcs->dAcc[0]    = 0;
	mcs->dAcc[1]    = 0;
	mcs->vAcc[0]    = 0;
	mcs->vAcc[1]    = 0;
	mcs->E          = 10E10;
	mcs->E6         = -5E10;
	mcs->E12        = 5E10;
	mcs->l          = mcs->N;
	mcs->Vir        = 10E10;
	mcs->Vir6       = -5E10;
	mcs->Vir12      = 5E10;
	mcs->md         = 0;

	mcs->cf = fopen("config.dat.mcs","w");
        mcs->tf = fopen("thermo.dat.mcs","w");
        fprintf(mcs->tf,"Step  Energy  Energy^2    l     l^2     Virial  Virial^2\n");
	fflush(mcs->tf);
        
        mcs->lA         = 0;
        mcs->lSA        = 0;
        mcs->EA         = 0;
        mcs->ESA        = 0;
        mcs->VirA       = 0;
        mcs->VirSA      = 0;

        mcs->dAcc[0]    = 0;
        mcs->dAcc[1]    = 0;
        mcs->vAcc[0]    = 0;
        mcs->vAcc[1]    = 0;

        mcs->numPairs   = ((mcs->N-1)*mcs->N)/2;

        mcs->r          = (double *) malloc(mcs->N*sizeof(double));
        mcs->rij        = (double *) malloc(mcs->numPairs*sizeof(double));
        rijTest         = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ij      = (double *) malloc(mcs->numPairs*sizeof(double));
        e12ijTest       = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ij       = (double *) malloc(mcs->numPairs*sizeof(double));
        e6ijTest        = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir12ij    = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir6ij     = (double *) malloc(mcs->numPairs*sizeof(double));

        mcs->rTrial         = (double *) malloc(mcs->N*sizeof(double));
        mcs->rijTrial       = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ijTrial      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir12ijTrial   = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir6ijTrial    = (double *) malloc(mcs->numPairs*sizeof(double));

        mcs->iii            = (unsigned long int *) malloc(mcs->numPairs*sizeof(double));
        mcs->jjj            = (unsigned long int *) malloc(mcs->numPairs*sizeof(double));
        mcs->indind         = (unsigned long int **) malloc((mcs->N-1)*sizeof(unsigned long int *));
        
	ind = 0;
        for (ii = 0; ii < mcs->N-1; ii++) {
                (*mcs).indind[ii] = (unsigned long int *) malloc((mcs->N-1-ii)*sizeof(unsigned long int));
                for (jj = ii + 1; jj < mcs->N; jj++) {
                        dj = jj - ii - 1;
                        mcs->indind[ii][dj] = ind;
                        ind++;
                }
        }

        mcs->slrho = -1;
        mcs->rhol = (int *) malloc(mcs->rhonb*sizeof(int));
        mcs->rhoA = (int *) malloc(mcs->rhonb*sizeof(int));
        mcs->rhoM = (double *) malloc(mcs->rhonb*sizeof(double));
        mcs->rhof = fopen("rho.dat.mcs","w");
        for (ii = 0; ii < mcs->N; ii++) {
                mcs->rhol[ii] = 0;
                mcs->rhoA[ii] = 0;
                mcs->rhoM[ii] = 0;
        }
        
        mcs->slg = -1;
        mcs->gl = (int **) malloc(mcs->gns*sizeof(int *));
        mcs->gA = (int **) malloc(mcs->gns*sizeof(int *));
        mcs->gM = (double **) malloc(mcs->gns*sizeof(double *));
        mcs->gf = (FILE **) malloc(mcs->gns*sizeof(FILE *));
        for (ii = 0; ii < mcs->gns; ii++) {
                mcs->gl[ii] = (int *) malloc(mcs->gnb*sizeof(int));
                mcs->gA[ii] = (int *) malloc(mcs->gnb*sizeof(int));
                mcs->gM[ii] = (double *) malloc(mcs->gnb*sizeof(double));
                sprintf(gfstr,"g%lu.dat.mcs",ii);
                mcs->gf[ii] = fopen(gfstr,"w");
                for (jj = 0; jj < mcs->gnb; jj++) {
                        mcs->gl[ii][jj] = 0;
                        mcs->gA[ii][jj] = 0;
                        mcs->gM[ii][jj] = 0;
                } 
        }       
        
        ind = 0;
        for (ii = 0; ii < mcs->N; ii++) {
                mcs->r[ii]        = -mcs->l/2 + (ii+0.5)*(mcs->l/mcs->N);
                for (jj = ii+1; jj < mcs->N; jj++) {
                        mcs->iii[ind] = ii;
                        mcs->jjj[ind] = jj;
                        ind++; 
                }
        }               
        for (ind = 0; ind < mcs->numPairs; ind++) {
                mcs->rij[ind] = mcs->r[mcs->jjj[ind]] - mcs->r[mcs->iii[ind]];
        }
        
        mcs_fgrho(mcs);
        mcs_ugrho(mcs);

	mcs->rangen = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(mcs->rangen,mcs->seed);
        fflush(stdout);
	
	return mcs;
}


void mcsPrintStep(struct MCState *mcs) {
	printf("Step: %lu...\n",mcs->sn);
}


void mcs_fad(struct MCState *mcs,unsigned long int *nm,double *d) {
	// Full attempt, displacement
	unsigned long int ii,jj,ind;
	double rij1,rij3,rij6,rij7,rij12,rij13,ran;
	
	

	#pragma omp single
	{
	if ( nm == NULL) {
		mcs->nm = gsl_rng_uniform_int(mcs->rangen,mcs->N);
	}
	else {
		mcs->nm = *nm;
	}
	if ( d == NULL ) {
		mcs->rn = gsl_rng_uniform(mcs->rangen);
	}
	else {
		mcs->rn = *d;
	}

	for (ind = 0; ind < mcs->N; ind++) {
		mcs->rTrial[ind] = mcs->r[ind];
	}
	for (ind = 0; ind < mcs->numPairs; ind++) {
		mcs->rijTrial[ind] = mcs->rij[ind];
	}

	mcs->md = (mcs->rn - 0.5)*2*mcs->maxStep;
	//fprintf(stdout,"Attempting to move particle %lu by %lf...",nm,md);
	mcs->rTrial[mcs->nm] = mcs->r[mcs->nm] + mcs->md;
	mcs->E6Trial = 0;
	mcs->E12Trial = 0;
	mcs->ETrial = 0;
	mcs->Vir6Trial = 0;
	mcs->Vir12Trial = 0;
	mcs->VirTrial = 0;
//	printf("Particle to move: %lu,  Move distance: %.5f...",nm,md);
	if (fabs(mcs->rTrial[mcs->nm]) > (mcs->l/2.0)) {
		mcs->bndchk = 1;
		//printf("Wall collision. Displacement rejected\n");
		mcs->dAcc[1]++;
	}
	else {
		mcs->bndchk = 0;
	}
	fflush(stdout);
	}
	
	mcsE6Trial     = mcs->E6Trial;
	mcsE12Trial    = mcs->E12Trial;
	mcsVir6Trial   = mcs->Vir6Trial;
	mcsVir12Trial  = mcs->Vir12Trial;
	
	if (mcs->bndchk == 0) {
		#pragma omp for private(ii,jj,ind,rij1,rij3,rij6,rij7,rij12,rij13) \
				reduction(+:mcsE6Trial,mcsE12Trial,mcsVir6Trial,mcsVir12Trial)
		for (ind = 0; ind < mcs->numPairs; ind++) {
			ii            = mcs->iii[ind];
			jj            = mcs->jjj[ind];
			mcs->rijTrial[ind] = mcs->rTrial[jj] - mcs->rTrial[ii];
//			printf("ii,jj,rij = %lu,%lu, %.5G\n",ii,jj,rijTrial[ind]);
			if (mcs->nbn >= 0 && jj-ii <= mcs->nbn) {
				mcs->e6ijTrial[ind]     = 0;
				mcs->e12ijTrial[ind]    = 0;
				mcs->vir6ijTrial[ind]   = 0;
				mcs->vir12ijTrial[ind]  = 0;	
			}
			else {
				rij1          = mcs->rijTrial[ind];
				rij3          = rij1*rij1*rij1;
				rij6          = 1/(rij3*rij3);
				rij7          = rij6/rij1;
				rij12         = rij6*rij6;
				rij13         = rij12/rij1;
	
				mcs->e6ijTrial[ind]     = 4*rij6;
				mcs->e12ijTrial[ind]    = 4*rij12;
				mcs->vir6ijTrial[ind]   = 24/mcs->l*rij6;
				mcs->vir12ijTrial[ind]  = 48/mcs->l*rij12;
			
				mcsE6Trial     += mcs->e6ijTrial[ind];
				mcsE12Trial    += mcs->e12ijTrial[ind];
				mcsVir6Trial   += mcs->vir6ijTrial[ind];
				mcsVir12Trial  += mcs->vir12ijTrial[ind];
			}

		}
//		printf("Done with for loop\n");
		
		mcs->E6Trial     = mcsE6Trial;
		mcs->E12Trial    = mcsE12Trial;
		mcs->Vir6Trial   = mcsVir6Trial;
		mcs->Vir12Trial  = mcsVir12Trial;
		
		
		#pragma omp single
		{
		mcs->ETrial = mcs->E12Trial - mcs->E6Trial;

		if (mcs->ETrial > mcs->E ) {
//			printf("ETrial: %lf...*Ep: %.5G...",ETrial,*Ep);
			ran = gsl_rng_uniform(mcs->rangen);
			mcs->bf = exp((mcs->E-mcs->ETrial)/mcs->T);
		}

		if (mcs->ETrial < mcs->E || mcs->bf > ran) {
			//printf("Displacement accepted...E = %.8G\n",ETrial);
			mcs->dAcc[0]++;
			mcs->E     = mcs->ETrial;
			mcs->E6    = mcs->E6Trial;
			mcs->E12   = mcs->E12Trial;
			mcs->Vir6  = mcs->Vir6Trial;
			mcs->Vir12 = mcs->Vir12Trial;
			mcs->Vir = mcs->N*mcs->T/mcs->l + mcs->Vir12 - mcs->Vir6;
			mcs->r[mcs->nm] = mcs->rTrial[mcs->nm];
			fflush(stdout);
			for (ind = 0; ind < mcs->numPairs; ind++) {
				mcs->rij[ind]     = mcs->rijTrial[ind];
				mcs->e6ij[ind]    = mcs->e6ijTrial[ind];
				mcs->e12ij[ind]   = mcs->e12ijTrial[ind];
				mcs->vir6ij[ind]  = mcs->vir6ijTrial[ind];
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
			}
		}
		else {
			//printf("Displacement rejected\n");
			mcs->dAcc[1]++;
		}
		//printf("Check 9...Thread no.: %d\n",omp_get_thread_num());
		}
	}
	//printf("Exiting fad...Thread no.: %d\n",omp_get_thread_num());
	fflush(stdout);
	#pragma omp barrier
}

int mcs_printCoords(struct MCState *mcs) {
        unsigned long int ind;
        fprintf(mcs->cf,"%lu\nStep no.: %lu  Box length: %.5f\n",mcs->N,mcs->sn,mcs->l);
        fflush(mcs->cf);
        for(ind = 0; ind < mcs->N; ind++) {
                fprintf(mcs->cf,"%lu  0.0  0.0  %.8G\n",ind+1,mcs->r[ind]);
        }
        fflush(mcs->cf);
	mcs->slcp = mcs->sn;
        return 0;
}


int mcs_printRho(struct MCState *mcs) {


        unsigned long int rb,ns;
        // if (fprintf(rhof,"") < 0 ) {
        // perror("Error printing to rhof.dat!\n");
        // }
        ns = (mcs->sn - mcs->slrho);
        //printf("printRho.  sn = %lu  rhof = %p  ns = %lu\n",sn,rhof,ns);
        //fflush(stdout);
        //printf("rhof new line\n");
        fprintf(mcs->rhof,"%lu",mcs->sn);
        //fflush(rhof);
        for (rb = 0; rb < mcs->rhonb; rb++) {
          //printf("## rhoM[%lu]: %p rhoA[%lu]: %p gA: %p %p gA[0]: %p %p  gA[%u]: %p %p\n",rb,&rhoM[rb],rb,&rhoA[rb],gA,&(gA[0]),gA[0],&(gA[0][0]),gns-1,gA[gns-1],&(gA[gns-1][0]));
          //fflush(stdout);
          //printf("rb = %lu\n",rb);
          //fflush(stdout);
          mcs->rhoM[rb] = (double) mcs->rhoA[rb]/ns/mcs->rbw;
          fprintf(mcs->rhof," %.8G",mcs->rhoM[rb]);
          //fflush(stdout);
          mcs->rhoA[rb] = 0;
        }
        fprintf(mcs->rhof,"\n");
        fflush(mcs->rhof);
        mcs->slrho = mcs->sn;

        return 0;
}

int mcs_printG(struct MCState *mcs) {
  
	int gs;
	unsigned long int gb,ns;
	
	#pragma omp for private(gs,gb,ns)
	for (gs = 0; gs < mcs->gns; gs++) {
		fprintf(mcs->gf[gs],"%lu",mcs->sn);
		ns = (mcs->sn - mcs->slg);

		for (gb = 0; gb < mcs->gnb; gb++) {
			mcs->gM[gs][gb] = (double) mcs->gA[gs][gb]/ns/mcs->gsw/mcs->gbw;
			fprintf(mcs->gf[gs]," %.8G",mcs->gM[gs][gb]);
			mcs->gA[gs][gb] = 0;
		}
		fprintf(mcs->gf[gs],"\n");
		fflush(mcs->gf[gs]);
	}
	mcs->slg = mcs->sn;

	return 0;

}
  
  



int mcs_fgrho(struct MCState *mcs) {

	unsigned long int ii,jj,ind;
	long int rb,gs1,gs2,gb;
	
	#pragma omp for
	for (rb = 0; rb < mcs->rhonb; rb++) {
		mcs->rhol[rb] = 0;
	}

	#pragma omp for
	for(ii = 0; ii < mcs->N; ii++) {
		//printf("ii = %lu\n",ii);
		rb = (long int) floor(mcs->r[ii]/mcs->rbw + mcs->rhonb/2.0);
		if (rb >= 0 && rb < mcs->rhonb) {
			#pragma omp atomic
			mcs->rhol[rb]++;
		}
	}
	
	//printf("Inside fgrho...Thread no.: %d...gns,gnb = %d, %lu\n",omp_get_thread_num(),gns,gnb);
	fflush(stdout);
	#pragma omp barrier
	// #pragma omp single
	// {
	
	#pragma omp for
	for (gs1 = 0; gs1 < mcs->gns; gs1++) {
		for (gb = 0; gb < mcs->gnb; gb++) {
			mcs->gl[gs1][gb] = 0;
		}
	}
	
	#pragma omp for private(ind,ii,jj,gs1,gs2,gb)
	//printf("gA: %p %p gA[0]: %p %p  gA[%u]: %p %p\n",gA,&(gA[0]),gA[0],&(gA[0][0]),gns-1,gA[gns-1],&(gA[gns-1][0]));
	for (ind = 0; ind < mcs->numPairs; ind++) {
		ii = mcs->iii[ind];
		jj = mcs->jjj[ind];
		gs1 = (long int) floor(mcs->r[ii]/mcs->gsw + mcs->gns/2.0);
		gs2 = (long int) floor(mcs->r[jj]/mcs->gsw + mcs->gns/2.0);
		gb = (long int) floor(mcs->rij[ind]/mcs->gbw);
		//printf("ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
		//fflush(stdout);
		if(gb < mcs->gnb) {
			if (gs1 >= 0 && gs1 < (long int) mcs->gns) {
				//printf("A ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
				//fflush(stdout);
				#pragma omp atomic
				mcs->gl[gs1][gb]++;
			}
			else {
				//printf("B ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
				//fflush(stdout);
				
			}
			
			if (gs2 >= 0 && gs2 < (long int) mcs->gns) {
				//printf("C ind = %lu  ii = %lu  jj = %lu  gs1 = %ld  gs2 = %ld  gb = %ld\n",ind,ii,jj,gs1,gs2,gb);
				//fflush(stdout);
				//printf("Accumulating to gA[%ld][%lu] at ",gs2,gb);
				//fflush(stdout);
				//printf("gA: %p  gA[%ld][%lu]: %p...",gA,gs2,gb,&(gA[gs2][gb]));
				//fflush(stdout);
				#pragma omp atomic
				mcs->gl[gs2][gb]++;
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

  return 0;
}



int mcs_ugrho(struct MCState *mcs) {

	unsigned long int rb,gb;
	unsigned int gs;
	
	#pragma omp for
	for (rb = 0; rb < mcs->rhonb; rb++) {
		mcs->rhoA[rb] += mcs->rhol[rb];
	}

	#pragma omp for private(gs,gb)
	for (gs = 0; gs < mcs->gns; gs++) {
		for (gb = 0; gb < mcs->gnb; gb++) {
			mcs->gA[gs][gb] += mcs->gl[gs][gb];
		}
	}

  return 0;
}



int mcs_qad(struct MCState *mcs) {

	unsigned long int ii,jj,dj,ind;
	double rij1,rij3,rij6,rij7,rij12,rij13;
	//double d;
	mcs->md = (mcs->rn - 0.5)*2*mcs->maxStep;
	
	#pragma omp single
	{
	mcs->rTrial[mcs->nm] = mcs->r[mcs->nm] + mcs->md;
	//fprintf(stdout,"Attempting to move particle %lu by %lf...",nm,md);
	}
	
	if (fabs(mcs->rTrial[mcs->nm]) > mcs->l/2.0) {
		#pragma omp single
		{
		mcs->dAcc[1]++;
		//printf("Wall collision. Move rejected...dRej = %lu",dAcc[1]);
		}
	}
	
	else {
		#pragma omp single
		{
		mcsdE6 = 0;
		mcsdE12 = 0;
		mcsdVir6 = 0;
		mcsdVir12 = 0;
		}
		
		//printf("nm = %lu...\n",nm);
		//fflush(stdout);	
		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) reduction(+:mcsdE6,mcsdE12,mcsdVir6,mcsdVir12) nowait
		for (ii = 0; ii < mcs->nm; ii++) {
		//	printf("ii = %lu\n",ii);
		//	fflush(stdout);	
			dj = mcs->nm - ii - 1;
			ind = mcs->indind[ii][dj];
			mcs->rijTrial[ind] = mcs->rij[ind] + mcs->md;
			if ( mcs->nbn >= 0 && dj >= mcs->nbn ) {
				mcs->e6ijTrial[ind] = 0;
				mcs->e12ijTrial[ind] = 0;
				mcs->vir6ijTrial[ind] = 0;
				mcs->vir12ijTrial[ind] = 0;
				//printf("nm,ii,e6,e6Trial,e12,e12Trial = %lu,%lu, %lf, %lf, %lf, %lf\n",nm,ii,e6ij[ind],e6ijTrial[ind],e12ij[ind],e12ijTrial[ind]);
			}
			else {
				rij1 = mcs->rijTrial[ind];
				rij3 = rij1*rij1*rij1;
				rij6 = 1/(rij3*rij3);
				rij7 = rij6/rij1;
				rij12 = rij6*rij6;
				rij13 = rij12/rij1;
				mcs->e6ijTrial[ind] = 4*rij6;
				mcs->e12ijTrial[ind] = 4*rij12;
				mcs->vir6ijTrial[ind] = 24/mcs->l*rij6;
				mcs->vir12ijTrial[ind] = 48/mcs->l*rij12;
				//printf("nm,ii,e6,e6Trial,e12,e12Trial = %lu,%lu, %lf, %lf, %lf, %lf\n",nm,ii,e6ij[ind],e6ijTrial[ind],e12ij[ind],e12ijTrial[ind]);
			}
			mcsdE6  = mcsdE6  - mcs->e6ij[ind]  + mcs->e6ijTrial[ind];
			mcsdE12 = mcsdE12 - mcs->e12ij[ind] + mcs->e12ijTrial[ind];
			mcsdVir6 = mcsdVir6 - mcs->vir6ij[ind] + mcs->vir6ijTrial[ind];
			mcsdVir12 = mcsdVir12 - mcs->vir12ij[ind] + mcs->vir12ijTrial[ind];
		}

		//#pragma omp single
		//{
		//	if (nm == 19) {
		//		fprintf(stdout,"(Intermediate) dE6, dE12 = %.8G, %.8G\n",dE6,dE12);
		//	}
		//}
		
		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) reduction(+:mcsdE6,mcsdE12,mcsdVir6,mcsdVir12)
		for (ii = mcs->nm + 1; ii < mcs->N; ii++) {
		//	printf("ii = %lu\n",ii);
		//	fflush(stdout);	
			dj = ii - mcs->nm - 1;
			ind = mcs->indind[mcs->nm][dj];
			mcs->rijTrial[ind] = mcs->rij[ind] - mcs->md;
			if (dj >= mcs->nbn) {
				mcs->e6ijTrial[ind]    = 0;
				mcs->e12ijTrial[ind]   = 0;
				mcs->vir6ijTrial[ind]  = 0;
				mcs->vir12ijTrial[ind] = 0;
				//printf("nm,ii,e6,e6Trial,e12,e12Trial = %lu,%lu, %lf, %lf, %lf, %lf\n",nm,ii,e6ij[ind],e6ijTrial[ind],e12ij[ind],e12ijTrial[ind]);
			}
			else {
				rij1 = mcs->rijTrial[ind];
				rij3 = rij1*rij1*rij1;
				rij6 = 1/(rij3*rij3);
				rij7 = rij6/rij1;
				rij12 = rij6*rij6;
				rij13 = rij12/rij1;
				mcs->e6ijTrial[ind] = 4*rij6;
				mcs->e12ijTrial[ind] = 4*rij12;
				mcs->vir6ijTrial[ind] = 24/mcs->l*rij6;
				mcs->vir12ijTrial[ind] = 48/mcs->l*rij12;
				//printf("nm,ii,e6,e6Trial,e12,e12Trial = %lu,%lu, %lf, %lf, %lf, %lf\n",nm,ii,e6ij[ind],e6ijTrial[ind],e12ij[ind],e12ijTrial[ind]);
			}
			mcsdE6  = mcsdE6  - mcs->e6ij[ind]  + mcs->e6ijTrial[ind];
			mcsdE12 = mcsdE12 - mcs->e12ij[ind] + mcs->e12ijTrial[ind];
			mcsdVir6 = mcsdVir6 - mcs->vir6ij[ind] + mcs->vir6ijTrial[ind];
			mcsdVir12 = mcsdVir12 - mcs->vir12ij[ind] + mcs->vir12ijTrial[ind];
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
		mcs->dE12 = mcsdE12;
		mcs->dE6 = mcsdE6;
		mcs->dE = mcs->dE12 - mcs->dE6;
		//printf("dE = %.8G...",dE);
		//fflush(stdout);	
		}
		
		if (mcs->dE > 0 ) {
		#pragma omp single
		{
			//printf("ETrial: %lf...*Ep: %.5G...",ETrial,*Ep);
			ran = gsl_rng_uniform(mcs->rangen);
			mcs->bf = exp(-mcs->dE/mcs->T);
		}
		}
		
		if (mcs->dE <= 0 || mcs->bf > ran) {
			#pragma omp single
			{
			mcs->dAcc[0]++;
                        //printf("Boltzmann factor: %lf...random: %lf...",bf,ran);
			//printf("Displacement accepted...dAcc = %lu...",dAcc[0]);
			mcs->E      += mcs->dE;
			mcs->E6     += mcs->dE6;
			mcs->E12    += mcs->dE12;
			mcs->Vir6   += mcs->dVir6;
			mcs->Vir12  += mcs->dVir12;
			mcs->Vir    += mcs->dVir12 - mcs->dVir6;
			mcs->r[mcs->nm] = mcs->rTrial[mcs->nm];
			}
			
			#pragma omp for private(jj,dj,ind) nowait
			for (jj = 0; jj < mcs->nm; jj++) {
				//printf("jj = %lu\n",jj);
				//fflush(stdout);	
				dj = mcs->nm - jj - 1;
				ind = mcs->indind[jj][dj];
				//printf("1. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->rij[ind] = mcs->rijTrial[ind];
				//printf("2. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->e6ij[ind] = mcs->e6ijTrial[ind];
				//printf("3. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->e12ij[ind] = mcs->e12ijTrial[ind];
				//printf("4. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
				//printf("5. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
				//printf("6. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				//fprintf(tf,"jj,nm,dj,rij = %lu,%lu,%lu,%lf\n",jj,nm,dj,rij[ind]);
			}
			#pragma omp for private(jj,dj,ind)
			for (jj = mcs->nm + 1; jj < mcs->N; jj++) {
				//printf("7. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				dj = jj - mcs->nm - 1;
				//printf("8. nm,jj,dj,ind = %lu,%lu,%lu,%lu\n",nm,jj,dj,ind);
				//fflush(stdout);	
				ind = mcs->indind[mcs->nm][dj];
				//printf("9. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->rij[ind] = mcs->rijTrial[ind];
				//printf("10. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->e6ij[ind] = mcs->e6ijTrial[ind];
				//printf("11. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->e12ij[ind] = mcs->e12ijTrial[ind];
				//printf("12. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
				//printf("13. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
				//printf("14. jj,dj,ind = %lu,%lu,%lu\n",jj,dj,ind);
				//fflush(stdout);	
				//fprintf(tf,"jj,nm,dj,rij = %lu,%lu,%lu,%lf\n",jj,nm,dj,rij[ind]);
				fflush(mcs->tf);
			}
			
			mcs_qagrho(mcs);
			
		/*#pragma omp single
		{
			//printf("Done accepting displacement\n");
			fflush(stdout);	
		}*/
		}
		else {
		#pragma omp single
		{
			mcs->dAcc[1]++;
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
	
	mcs_ugrho(mcs);
	//printf("Exiting fad...Thread no.: %d\n",omp_get_thread_num());
	fflush(stdout);
	#pragma omp barrier


	return 0;
}



int mcs_qav(struct MCState *mcs) {
	
	unsigned long int ii,ind;
	double dQ,ran,lRat1,lRat3,lRat6,lRat7,lRat12,lRat13;
	
	#pragma omp single
	{
	mcs->dl = (mcs->rn-0.5)*2*mcs->maxdl;
	//fprintf(stdout,"Attempting to change volume by %lf...",dl);
	//printf("E6 = %.5G...E12 = %.5G...E = %.5G...",E6,E12,E);
	lRat1 = (mcs->l+mcs->dl)/mcs->l;
	lRat3 = lRat1*lRat1*lRat1;
	lRat6 = 1/(lRat3*lRat3);
	lRat12 = lRat6*lRat6;
	mcs->E6Trial = lRat6*mcs->E6;
	mcs->E12Trial = lRat12*mcs->E12;
	mcs->dE = mcs->E12Trial - mcs->E6Trial - mcs->E;
	//mcs->bf = exp(-(dE + P*dl)/T)*pow(lRat1,N);
	mcs->bf   = exp(-(mcs->dE + mcs->P*mcs->dl)/mcs->T + mcs->N*log(lRat1));
	//printf("E6Trial = %.5G...E12Trial = %.5G...dE = %.5G...Boltzmann factor: %lf",E6Trial,E12Trial,dE,bf);
	//printf("dE = %.5G...Boltzmann factor: %lf...",dE,bf);

	if (mcs->bf < 1.0 ) {
//		printf("ETrial: %lf...*Ep: %.5G...",ETrial,*Ep);
		ran = gsl_rng_uniform(mcs->rangen);
		//printf("dE: %lf...Boltzmann factor: %lf...",dE,bf);
	}
	}
	
	if (mcs->bf >= 1.0 || mcs->bf > ran) {
		#pragma omp single
		{
		//printf("random: %lf...",ran);
		mcs->vAcc[0]++;
		//printf("Volume change accepted...vAcc = %lu",vAcc[0]);
		mcs->E6  = mcs->E6Trial;
		mcs->E12 = mcs->E12Trial;
		mcs->E      = mcs->E12 - mcs->E6;
		
		mcs->l = mcs->l + mcs->dl;
		
		lRat7 = lRat6/lRat1;
		lRat13 = lRat12/lRat1;
		mcs->Vir6 = lRat7*mcs->Vir6;
		mcs->Vir12 = lRat13*mcs->Vir12;
		mcs->Vir    = mcs->N*mcs->T/mcs->l + mcs->Vir12 - mcs->Vir6;
		}
		
		#pragma omp single private(ii) nowait
		for (ii = 0; ii < mcs->N; ii++) {
			mcs->r[ii] = lRat1*mcs->r[ii];
		}
		
		#pragma omp for private(ind)
		for (ind = 0; ind < mcs->numPairs; ind++) { 
			mcs->rij[ind]     = lRat1*mcs->rij[ind];
			mcs->e6ij[ind]    = lRat6*mcs->e6ij[ind];
			mcs->e12ij[ind]   = lRat12*mcs->e12ij[ind];
			mcs->vir6ij[ind]  = lRat7*mcs->vir6ij[ind];
			mcs->vir12ij[ind] = lRat13*mcs->vir12ij[ind];
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
		mcs_fgrho(mcs);
	}
	
	else {
		#pragma omp single
		{
		//printf("random: %lf...",ran);
		mcs->vAcc[1]++;
		//printf("Volume change rejected...vRej = %lu",vAcc[1]);
		}
	}
	mcs_ugrho(mcs);
	fflush(stdout);

	return 0;
}






unsigned long int mcs_incrementStep(struct MCState *mcs) {
  int flag = 0;
  const int stepPrintInterval = 1;
  
  if (mcs->sn == mcs->numSteps) {
    return 0;
  }
  else {
    #pragma omp barrier
    #pragma omp single
    {
    mcs->sn++;
    
    if (mcs->sn % stepPrintInterval == 0) {
      printf("Step: %lu...\n",mcs->sn);
      //fflush(stdout);
    }
    }
    return mcs->sn;
  }
  #pragma omp barrier
}







int mcs_Step(struct MCState *mcs) {
  #pragma omp barrier
  #pragma omp single
  {
  mcs->nm = gsl_rng_uniform_int(mcs->rangen,mcs->N+1);
  mcs->rn = gsl_rng_uniform(mcs->rangen);
  printf("nm: %lu, rn: %.5G...\n",mcs->nm,mcs->rn);
  fflush(stdout);
  }
  #pragma omp barrier

  if (mcs->nm < mcs->N) {
    printf("Executing mcs_qad(mcs)...\n");
    fflush(stdout);
    mcs_qad(mcs);
  }
  else {
    printf("Executing mcs_qav(mcs)...\n");
    fflush(stdout);
    mcs_qav(mcs);
  }
  
  if (mcs_isECheck(mcs) ) {
    mcs_ECheck(mcs);
  }
  
  mcs_updateThermo(mcs);

  return 0;

}


int mcs_isCoordPrint(struct MCState *mcs) {
  if ( mcs->sn % mcs->cpi == 0 ) {
    return 1;
  }
  else {
    return 0;
  }
}


int mcs_isThermoPrint(struct MCState *mcs ) {
  if ( mcs->sn % mcs->tpi == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}

int mcs_isRhoPrint(struct MCState *mcs) {
  if ( mcs->sn % mcs->rhopi == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}

int mcs_isGPrint(struct MCState *mcs) {
  if ( mcs->sn % mcs->gpi == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}

int mcs_isECheck(struct MCState *mcs) {
  if ( mcs->sn % mcs->eci == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}

int mcs_isMaxDisAdjust(struct MCState *mcs) {
  if ( mcs->sn % mcs->mdai == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}

int mcs_isMaxDVAdjust(struct MCState *mcs) {
  if ( mcs->sn % mcs->mvai == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}


int mcs_printThermo(struct MCState *mcs) {
	double lM,lSM,EM,ESM,VirM,VirSM;
	unsigned long int ss;

	ss = mcs->sn - mcs->sltp;

	lM       = mcs->lA/ss;
	lSM      = mcs->lSA/ss;
	EM       = mcs->EA/ss;
	ESM      = mcs->ESA/ss;
	VirM     = mcs->VirA/ss;
	VirSM    = mcs->VirSA/ss;
//	printf("Check 3...sn: %lu...EM: %.5G...ESM: %.5G...lM: %.5G...lSM: %.5G\n",sn,EM,ESM,lM,lSM);
	
	fprintf(mcs->tf,"%lu\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\n",mcs->sn,EM,ESM,lM,lSM,VirM,VirSM);
	fflush(mcs->tf);
	mcs->lA     = 0;
	mcs->lSA    = 0;
	mcs->EA     = 0;
	mcs->ESA    = 0;
	mcs->VirA   = 0;
	mcs->VirSA  = 0;
	mcs->sltp   = mcs->sn;

	return 0;
}


int mcs_updateThermo(struct MCState *mcs) {

    mcs->lA = mcs->lA + mcs->l;
    mcs->lSA = mcs->lSA + mcs->l*mcs->l; 
    mcs->EA = mcs->EA + mcs->E; 
    mcs->ESA = mcs->ESA + mcs->E*mcs->E; 
    mcs->VirA = mcs->VirA + mcs->Vir; 
    mcs->VirSA = mcs->VirSA + mcs->Vir*mcs->Vir;

    return 0;

}



int mcs_ECheck(struct MCState *mcs) {
  int resultFlag;
  unsigned long int iq,ii,jj;
  double rij1,rij3,rij6,rij12;
  
  #pragma omp for private(iq,ii,jj,rij1,rij3,rij6,rij12) reduction(+:mcsETest)
  for (iq = 0; iq < mcs->numPairs; iq++) {
    ii = mcs->iii[iq];
    jj = mcs->jjj[iq];
    
    rijTest[iq] = mcs->r[jj] - mcs->r[ii];
    if (mcs->nbn >= 0 && jj-ii > mcs->nbn) {
      e12ijTest[iq] = 0;
      e6ijTest[iq]  = 0;
    }
    else {
      rij1 = rijTest[iq];
      rij3 = rij1*rij1*rij1;
      rij6 = 1/(rij3*rij3);
      rij12 = rij6*rij6;
      e12ijTest[iq] = 4*rij12;
      e6ijTest[iq] = 4*rij6;
      mcsETest += e12ijTest[iq] - e6ijTest[iq];
    }
  }
  #pragma omp single
  {
  if (fabs(mcsETest - mcs->E) > 0.001) {
    resultFlag = 1;
    fprintf(stdout,"\n\n\n### !!! $$$ Energy discrepancy!!! ######### !!!!!!!!!!! $$$$$"
         "$$$$$ !!!!!!!!\nE = %.8G, ETest = %.8G, Resetting energy "
         "to that of full calculation (E = %.8G)\n\n",mcs->E,mcsETest,mcsETest);
    mcs->E = mcsETest;
  }
  else {
    resultFlag = 0;
    fprintf(stdout,"####### Energy was just verified ########  E = %.8G ######## ETest = %.8G #######\n",mcs->E,mcsETest);
  }
  }
  return resultFlag;
}
  

int mcs_maxDisAdjust(struct MCState *mcs) {
  static unsigned long int dAErrNtot = 0;
  static double idealRatio = 0.5;
  double actualRatio;
  
  if ( (mcs->dAcc[0] + mcs->dAcc[1] - dAErrNtot) != 0 && (mcs->dAcc[0] + mcs->dAcc[1]) % mcs->mdai == 0 ) {
    dAErrNtot = mcs->dAcc[0] + mcs->dAcc[1];
    actualRatio = (double) mcs->dAcc[0] / (mcs->dAcc[0] + mcs->dAcc[1]);
    mcs->maxStep = mcs->maxStep*log(0.672924*idealRatio + 0.0644284)/log(0.672924*(actualRatio + 0.0644284));
    printf("Step: %lu  dAcc: %lu,%lu   maxStep: %.5G\n",mcs->sn,mcs->dAcc[0],mcs->dAcc[1],mcs->maxStep);
  }
    
  return 0;
}

int mcs_maxDVAdjust(struct MCState *mcs) {
  static unsigned long int vAErrNtot = 0;
  static double idealRatio = 0.5;
  double actualRatio;
  
  if ( (mcs->vAcc[0] + mcs->vAcc[1] - vAErrNtot) != 0 && (mcs->vAcc[0] + mcs->vAcc[1]) % mcs->mvai == 0 ) {
    vAErrNtot = mcs->vAcc[0] + mcs->vAcc[1];
    actualRatio = (double) mcs->vAcc[0] / (mcs->vAcc[0] + mcs->vAcc[1]);
    mcs->maxdl = mcs->maxdl*log(0.672924*idealRatio + 0.0644284)/log(0.672924*(actualRatio + 0.0644284));
    printf("Step: %lu  vAcc: %lu,%lu   maxStep: %.5G\n",mcs->sn,mcs->vAcc[0],mcs->vAcc[1],mcs->maxdl);
  }
  
  return 0;
}


int mcs_printE(struct MCState *mcs) {
  printf("\nE = %.8G\n",mcs->E);
  return 0;
}

int mcs_printAcc(struct MCState *mcs) {
  printf("Accepted/Rejected: Displacements VolumeChanges\n              \
          %lu/%lu          %lu/%lu\n",mcs->dAcc[0],mcs->dAcc[1],mcs->vAcc[0], \
          mcs->vAcc[1]);
  return 0;
}






//void fav(dVECTOR r, dVECTOR rTrial, dVECTOR rij, dVECTOR rijTrial, dVECTOR e6ij, dVECTOR e6ijTrial, dVECTOR e12ij, dVECTOR e12ijTrial,
//		dVECTOR vir6ij, dVECTOR vir6ijTrial, dVECTOR vir12ij, dVECTOR vir12ijTrial, double *lp,double *Ep, double *Virp, double rn) {
int mcs_fav(struct MCState *mcs) {
	unsigned long int ii,jj,ind;
	double lRat1,ran;
	double rij1a,rij3a,rij6a,rij7a,rij12a,rij13a;

	#pragma omp single
	{
	mcs->bndchk = 0;
	mcsETrial = 0;
	mcsVirTrial = 0;
	mcs->dl = (mcs->rn-0.5)*2*mcs->maxdl;
//	printf("Volume change attempt. Linitial = %.5G, dL = %.5G...",*lp,dl);
	lRat1 = (mcs->l + mcs->dl) / mcs->l;
//	printf("dl = %.5G, lRat1 = %.5G\n",dl,lRat1);
	for (ind = 0; ind < mcs->N; ind++) {
		mcs->rTrial[ind] = mcs->r[ind]*lRat1;
	}

	}

	#pragma omp for private(ii,jj,ind,rij1a,rij3a,rij6a,rij7a,rij12a,rij13a) reduction (+:mcsETrial,mcsVirTrial)
	for (ind = 0; ind < mcs->numPairs; ind++) {
		ii            = mcs->iii[ind];
		jj            = mcs->jjj[ind];

//		printf("ii,jj,rij,= %lu,%lu, %.5G\n",ii,jj,rij[ind]);

		rij1a          = mcs->rTrial[jj] - mcs->rTrial[ii];
		if ( mcs->nbn >= 0 && jj-ii > mcs->nbn) {
//			printf("ii,jj,rij,= %lu,%lu, %.5G\n",ii,jj,rij[ind]);
			mcs->e6ijTrial[ind]       = 0;
			mcs->e12ijTrial[ind]      = 0;
			mcs->vir6ijTrial[ind]     = 0;
			mcs->vir12ijTrial[ind]    = 0;
//			printf("ii,jj,rij1,rij6,e6 = %lu,%lu, %.5G, %.5G, %.5G\n",ii,jj,rij[ind],rij6a,e6ij[ind]);
		}
		else {
//			printf("ii,jj,rij,= %lu,%lu, %.5G\n",ii,jj,rij[ind]);
			rij3a          = rij1a*rij1a*rij1a;
			rij6a          = 1/(rij3a*rij3a);
			rij7a          = rij6a/rij1a;
			rij12a         = rij6a*rij6a;
			rij13a         = rij12a/rij1a;
	
			mcs->e6ijTrial[ind]     = 4*rij6a;
			mcs->e12ijTrial[ind]    = 4*rij12a;
			mcs->vir6ijTrial[ind]     = 24/mcs->l*rij6a;
			mcs->vir12ijTrial[ind]    = 48/mcs->l*rij12a;
//			printf("ii,jj,rij1,rij6,e6 = %lu,%lu, %.5G, %.5G, %.5G\n",ii,jj,rij[ind],rij6a,e6ij[ind]);
	
			mcsETrial            += mcs->e12ijTrial[ind] - mcs->e6ijTrial[ind];
			mcsVirTrial          += mcs->vir12ijTrial[ind] - mcs->vir6ijTrial[ind];
		}
	}
			mcs->ETrial           = mcsETrial;
			mcs->VirTrial         = mcsVirTrial;

	#pragma omp single
	{

//	printf("ETrial: %lf...",ETrial);
	mcs->bf = exp(-(mcs->ETrial - mcs->E + mcs->P*mcs->dl) / mcs->T + mcs->N*log(lRat1));

	if (mcs->bf < 1.0) {
		ran = gsl_rng_uniform(mcs->rangen);
	}

	if (mcs->bf >= 1.0 || mcs->bf > ran) {
//		printf("Volume change accepted...");
		mcs->vAcc[0]++;
		mcs->l = mcs->l + mcs->dl;
		mcs->E = mcs->ETrial;
		mcs->Vir = mcs->VirTrial;
		for (ind = 0; ind < mcs->N; ind++) {
			mcs->r[ind] = mcs->rTrial[ind];
		}
		for (ind = 0; ind < mcs->numPairs; ind++) {
			mcs->rij[ind] = lRat1*mcs->rijTrial[ind];
			mcs->e6ij[ind] = mcs->e6ijTrial[ind];
			mcs->e12ij[ind] = mcs->e12ijTrial[ind];
			mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
			mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
		}
//		printf("Volume now: %.5G\n",*lp);
	}
	else {
		mcs->vAcc[1]++;
//		printf("Volume change rejected\n");
	}
	fflush(stdout);
	}
	
	return 0;
	#pragma omp barrier
}




int mcs_qagrho(struct MCState *mcs) {
	
	unsigned long int ii,jj;
	long int rb,gs,gb,gb1,gb2,gs2;
	int aa,bb;
		
	#pragma omp single
	{
	static long int rbn1,rbn2;
	rbn1 = (long int) floor((mcs->r[mcs->nm] - mcs->md)/mcs->rbw + mcs->rhonb/2.0);
	rbn2 = (long int) floor(mcs->r[mcs->nm]/mcs->rbw + mcs->rhonb/2.0);
	if (rbn1 >= 0 && rbn1 < mcs->rhonb) {
		mcs->rhol[rbn1]--;
	}
	if (rbn2 >= 0 && rbn2 < mcs->rhonb) {
		mcs->rhol[rbn2]++;
	}
	}
	
	#pragma omp single
	{
	gs11 = (long int) floor((mcs->r[mcs->nm] - mcs->md)/mcs->gsw + mcs->gns/2.0);
	gs12 = (long int) floor( mcs->r[mcs->nm]           /mcs->gsw + mcs->gns/2.0);
	}
	
	aa = 0;
	bb = 0;
	if (gs11 >= 0 && gs11 < mcs->gns) aa += 1;
	if (gs12 >= 0 && gs12 < mcs->gns) aa += 2;

	#pragma omp for private(ii,jj,gs2,gb1,gb2) nowait
	for (ii = 0; ii < mcs->nm; ii++) {
        	bb = aa;
	        gs2 = (long int) floor(mcs->r[ii]/mcs->gsw + mcs->gns/2.0);
		if (gs2 >= 0 && gs2 < mcs->gns) bb += 4;
		gb1 = (unsigned long int) floor(fabs((mcs->r[mcs->nm] - mcs->md - mcs->r[ii]))/mcs->gbw);
		if (gb1 >= 0 && gb1 < mcs->gnb) bb += 8;
		gb2 = (unsigned long int) floor(fabs((mcs->r[mcs->nm]      - mcs->r[ii]))/mcs->gbw);
		if (gb2 >= 0 && gb2 < mcs->gnb) bb += 16;
		
		if (bb % 2 >= 1 && bb % 16 >= 8) {
			#pragma omp atomic
			mcs->gl[gs11][gb1]--;
		}
		if (bb % 4 >= 2 && bb % 32 >= 16) {
			#pragma omp atomic
			mcs->gl[gs12][gb2]++;
		}
		if (bb % 8 >= 4 && bb % 16 >= 8) {
			#pragma omp atomic
			mcs->gl[gs2][gb1]--;
		}
		if (bb % 8 >= 4 && bb % 32 >= 16) {
			#pragma omp atomic
			mcs->gl[gs2][gb2]++;
		}
	}
	
	#pragma omp for private(ii,jj,gs2,gb1,gb2)
	for (ii = mcs->nm + 1; ii < mcs->N; ii++) {
        	bb = aa;
	        gs2 = (long int) floor(mcs->r[ii]/mcs->gsw + mcs->gns/2.0);
		if (gs2 >= 0 && gs2 < mcs->gns) bb += 4;
		gb1 = (unsigned long int) floor(fabs((mcs->r[mcs->nm] - mcs->md - mcs->r[ii]))/mcs->gbw);
		if (gb1 >= 0 && gb1 < mcs->gnb) bb += 8;
		gb2 = (unsigned long int) floor(fabs((mcs->r[mcs->nm]      - mcs->r[ii]))/mcs->gbw);
		if (gb2 >= 0 && gb2 < mcs->gnb) bb += 16;

		if (bb % 2 >= 1 && bb % 16 >= 8) {
			#pragma omp atomic
			mcs->gl[gs11][gb1]--;
		}
		if (bb % 4 >= 2 && bb % 32 >= 16) {
			#pragma omp atomic
			mcs->gl[gs12][gb2]++;
		}
		if (bb % 8 >= 4 && bb % 16 >= 8) {
			#pragma omp atomic
			mcs->gl[gs2][gb1]--;
		}
		if (bb % 8 >= 4 && bb % 32 >= 16) {
			#pragma omp atomic
			mcs->gl[gs2][gb2]++;
		}
	}
	
  return 0;
}






