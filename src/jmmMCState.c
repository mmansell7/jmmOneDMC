#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "jmmMCState.h"

double mcsE6Trial, mcsE12Trial, mcsVir6Trial, mcsVir12Trial;

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
  unsigned long int sn;            // Step number
  unsigned long int numSteps;      // Number of steps (total steps requested)
  unsigned long int cpi;           // Configuration print interval
  unsigned long int tpi;           // Thermo print interval
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
  gsl_rng * rangen;                  // Random number generator
};

int printMCP(struct MCState *mcs1) {
	printf("Printing Monte Carlo parameters...\n");
        printf("N: %lu\nP: %.5G\nT: %.5G\nnumSteps: %lu\nmaxStep: %.5G\nmax vol change: %.5G\n",mcs1->N,mcs1->P,mcs1->T,mcs1->numSteps,mcs1->maxStep,mcs1->maxdl);
        printf("Configuration print interval: %lu\nThermo print interval: %lu\nDensity bin width: %.5G\n",mcs1->cpi,mcs1->tpi,mcs1->rbw);
        printf("Number of density bins: %lu\nDensity print interval: %lu\ng(x) (or two-particle density) segment width: %.5G\n",mcs1->rhonb,mcs1->rhopi,mcs1->gsw);
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
	mcs->sn         = 0;
        mcs->numSteps   = inp.ns;
        mcs->cpi        = inp.cpi;
        mcs->tpi        = inp.tpi;
	mcs->slcp       = 0;
	mcs->sltp       = 0;
	mcs->numPairs   = (unsigned long int) ((double) (mcs->N-1)/2*mcs->N);
	mcs->nm         = 0;
	mcs->gpi        = inp.gpi;
	mcs->rhopi      = inp.rhopi;
	mcs->gnb        = inp.gnb;
	mcs->rhonb      = inp.rhonb;
	mcs->maxStep    = inp.maxStep;
	mcs->maxdl      = inp.maxdl;
	mcs->phi        = inp.phi;
	mcs->slg        = 0;
	mcs->slrho      = 0;
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
	mcs->lA         = 0;
	mcs->lSA        = 0;
	mcs->EA         = 0;
	mcs->ESA        = 0;
	mcs->VirA       = 0;
	mcs->VirSA      = 0;


	mcs->cf = fopen("config.dat","w");
        mcs->tf = fopen("thermo.dat","w");
        fprintf(mcs->tf,"Step  Energy  Energy^2    l     l^2     Virial  Virial^2\n");
	fflush(mcs->tf);
        
	mcs->slcp       = 0;
        mcs->sltp       = 0;
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
        mcs->E          = 10E10;
        mcs->l          = mcs->N;
        mcs->Vir        = 10E10;

        mcs->r          = (double *) malloc(mcs->N*sizeof(double));
        mcs->rij        = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ij      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ij       = (double *) malloc(mcs->numPairs*sizeof(double));
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

        mcs->slrho = 0;
        mcs->rhol = (int *) malloc(mcs->rhonb*sizeof(int));
        mcs->rhoA = (int *) malloc(mcs->rhonb*sizeof(int));
        mcs->rhoM = (double *) malloc(mcs->rhonb*sizeof(double));
        mcs->rhof = fopen("rho.dat","w");
        for (ii = 0; ii < mcs->N; ii++) {
                mcs->rhol[ii] = 0;
                mcs->rhoA[ii] = 0;
                mcs->rhoM[ii] = 0;
        }
        
        mcs->slg = 0;
        mcs->gl = (int **) malloc(mcs->gns*sizeof(int *));
        mcs->gA = (int **) malloc(mcs->gns*sizeof(int *));
        mcs->gM = (double **) malloc(mcs->gns*sizeof(double *));
        mcs->gf = (FILE **) malloc(mcs->gns*sizeof(FILE *));
        for (ii = 0; ii < mcs->gns; ii++) {
                mcs->gl[ii] = (int *) malloc(mcs->gnb*sizeof(int));
                mcs->gA[ii] = (int *) malloc(mcs->gnb*sizeof(int));
                mcs->gM[ii] = (double *) malloc(mcs->gnb*sizeof(double));
                sprintf(gfstr,"g%lu.dat",ii);
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
        
        //fgrho();
        //ugrho();

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
//			printf("ind,ii,jj,rijTrial,e6Trial,e12Trial,ETrial = %lu,%lu,%lu, %.5G, %.5G, %.5G, %.5G\n",ind,ii,jj,rijTrial[ind],e6ijTrial[ind],e12ijTrial[ind],ETrial);

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

