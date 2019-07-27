#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include "jmmMCState.h"
#include "pot.h"
#include "stdbool.h"

// Some file-global variables (available to all functions defined in this
//  file, but not to functions defined in other files) are necessary
//  to serve as reduction variables for OpenMP for loops. Local variables
//  cannot be used for this purpose.
static double mcsETest,mcsETrial,mcsE6Trial,mcsE12Trial,mcsVirTrial,mcsVir6Trial,mcsVir12Trial;
static double mcsdE6,mcsdE12,mcsdE,mcsdVir6,mcsdVir12,mcsdVir,ran;
static int resultFlag;
static double EUp,EDown,h,dlEstimate,lTryMin,lTryMax,b;
static double *rijTest,*eijTest,*e6ijTest,*e12ijTest;
static long int gs11,gs12;
double lRat1,lRat3,lRat6,lRat7,lRat12,lRat13;
static double params[10];

// IMPORTANT. Define the MCState struct. This struct holds the state of the system
//  and is what other functions operate on in order to read or modify that state.
//  Note that, by defining the struct in jmmMCState.c instead of jmmMCState.h, its
//  internal variables cannot be accessed except through the functions defined in 
//  this file. This provides a level of safety.
struct MCState {
  
  // Note: Use of typedefs to hide pointers is discouraged.  I have gone away from
  //  that practice.

  bool isRestart;                  // Boolean flag determining whether this is a restart run
                                   //   false: not restart
  
  unsigned long int N;             // Number of particles
  int nbn;                         // Number of neighbors with which each particle can
                                   //   interact.  -1 for no limit.
  unsigned long int sn;            // Step number
  unsigned long int numSteps;      // Number of steps (total steps requested)
  int relaxFlag;                   // 0: do not relax before run.  1: do relax
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
  double rhoSysA;                  // Accumulator for system overall (mean) density
  double rhoSysSA;                 // rhoSys^2 accumulator
  double lA;                       // l accumulator
  double lSA;                      // l^2 accumulator
  double EA;                       // E accumulator
  double ESA;                      // E^2 accumulator
  double lEA;                      // l*E accumulator
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
  double *eij;                   // Vector of pair contributions to E12
  double *eijTrial;              // e12ij for current trial
  double *e12ij;                   // Vector of pair contributions to E12
  double *e12ijTrial;              // e12ij for current trial
  double *e6ij;                    // Vector of pair contributions to E6
  double *e6ijTrial;               // e6ij for current trial 
  double *virij;                  // Vector of pair contributions to Vir6
  double *virijTrial;             // vir6ij for current trial
  double *vir12ij;                 // Vector of pair contributions to Vir12
  double *vir12ijTrial;            // vir12ij for current trial
  double *vir6ij;                  // Vector of pair contributions to Vir6
  double *vir6ijTrial;             // vir6ij for current trial

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
  
  
  char potStr[20];                 // String defining the potential type
  double potCutOff;
  void (*phi)(double *rij,double cutOff,         // Pointer to the
      void *params, double phi[6]);//   interparticle potential function
  int (*qad)(struct MCState *,     // Pointer a potential-specific, quicker,
     unsigned long int *nm,double *d);  // displacement trial function
  int (*qav)(struct MCState *);    // Pointer to a potential-specific, quicker,
                                   //   volume change trial function
};


// Print to stdout the (not necessarily) constant "parameters" defining the 
//  requested simulation conditions.
int printMCP(struct MCState *mcs1) {
        if (mcs1->isRestart) {
            printf("*** THIS IS A RESTART RUN ***\n");
        }
	printf("Printing Monte Carlo parameters...\n");
        printf("N: %lu\nP: %.5G\nT: %.5G\n",mcs1->N,mcs1->P,mcs1->T);
	printf("Potential: %s\n",mcs1->potStr);
        printf("Potential cut-off: %.5G\n",mcs1->potCutOff);
        if ( mcs1->nbn > 0 ) {
		printf("Number of neighbors with which each particle can interact: %u\n",mcs1->nbn);
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
        printf("Adjust max. displacement every %lu steps.\n",mcs1->mdai);
	printf("Adjust max. volume change every %lu steps.\n",mcs1->mvai);
	fflush(stdout);
	
	return 0;
}	


// Take an MCInput struct and use it to create and return a new MCState struct.
struct MCState * setupMCS(struct MCInput inp) {
    struct MCState *mcs = (struct MCState *) malloc(sizeof(struct MCState));

    unsigned long int ii,jj,dj,ind;
    char gfstr[20];
    
    // Initialize parameters provided by input struct.
    mcs->isRestart  = inp.isRestart;
    mcs->N          = inp.N;
    mcs->P          = inp.P;
    mcs->T          = inp.T;
    mcs->nbn        = inp.nbn;
    mcs->numSteps   = inp.ns;
    mcs->relaxFlag  = inp.relaxFlag;
    if ( mcs->isRestart && mcs->relaxFlag ) {
        printf("WARNING: running restart job with volume relaxation. Make sure you\n"
               "know what the heck you are doing!\n\n");
    }
    mcs->cpi        = inp.cpi;
    mcs->tpi        = inp.tpi;
    mcs->eci        = inp.eci;
    mcs->mdai       = inp.mdai;
    mcs->mvai       = inp.mvai;
    mcs->numPairs   = (unsigned long int) ((double) (mcs->N-1)/2*mcs->N);
    mcs->nm         = 0;
    mcs->gpi        = inp.gpi;
    mcs->rhopi      = inp.rhopi;
    mcs->gnb        = inp.gnb;
    mcs->rhonb      = inp.rhonb;
    mcs->maxStep    = inp.maxStep;
    mcs->maxdl      = inp.maxdl;
    mcs->rbw        = inp.rbw;
    mcs->gsw        = inp.gsw;
    mcs->gbw        = inp.gbw;
    mcs->gns        = inp.gns;
    mcs->seed       = inp.seed;
    mcs->potCutOff  = inp.potCutOff;
    strncpy(mcs->potStr,inp.potStr,20);
    if ( strncmp(mcs->potStr,"LJ",10) == 0 ) {
        mcs->phi        = &phiLJ;
        mcs->qad        = &qad2;
        if ( mcs->nbn < 0 ) {
            mcs->qav    = &qavLJ;
        }
        else {
            mcs->qav    = &fav;
        }
        mcs->E6         = -5E10;
        mcs->E12        = 5E10;
        mcs->Vir6       = -5E10;
        mcs->Vir12      = 5E10;
        mcs->eij        = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ij      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ij       = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->virij      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir12ij    = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir6ij     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->eijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ijTrial      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->virijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir12ijTrial   = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir6ijTrial    = (double *) malloc(mcs->numPairs*sizeof(double));
        eijTest         = (double *) malloc(mcs->numPairs*sizeof(double));
        e12ijTest       = (double *) malloc(mcs->numPairs*sizeof(double));
        e6ijTest        = (double *) malloc(mcs->numPairs*sizeof(double));
    }
    else if ( strncmp(mcs->potStr,"LJcut",10) == 0 ) {
        mcs->phi        = &phiLJcut;
        mcs->qad        = &qad2;
        mcs->qav        = &fav;
        mcs->E6         = -5E10;
        mcs->E12        = 5E10;
        mcs->Vir6       = -5E10;
        mcs->Vir12      = 5E10;
        mcs->eij        = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ij      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ij       = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->virij      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir12ij    = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir6ij     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->eijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ijTrial      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->virijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir12ijTrial   = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir6ijTrial    = (double *) malloc(mcs->numPairs*sizeof(double));
        eijTest         = (double *) malloc(mcs->numPairs*sizeof(double));
        e12ijTest       = (double *) malloc(mcs->numPairs*sizeof(double));
        e6ijTest        = (double *) malloc(mcs->numPairs*sizeof(double));
    }
    else if ( strncmp(mcs->potStr,"HARMONIC",10) == 0 ) {
        mcs->phi        = phiHarmonic;
        mcs->qad        = &qad2;
        mcs->qav        = &fav;
        mcs->E6         = -5E10;
        mcs->E12        = 5E10;
        mcs->Vir6       = -5E10;
        mcs->Vir12      = 5E10;
        mcs->eij        = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ij      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ij       = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->virij      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir12ij    = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir6ij     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->eijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e12ijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->e6ijTrial      = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->virijTrial     = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir12ijTrial   = (double *) malloc(mcs->numPairs*sizeof(double));
        mcs->vir6ijTrial    = (double *) malloc(mcs->numPairs*sizeof(double));
        eijTest         = (double *) malloc(mcs->numPairs*sizeof(double));
        e12ijTest       = (double *) malloc(mcs->numPairs*sizeof(double));
        e6ijTest        = (double *) malloc(mcs->numPairs*sizeof(double));
    }
    else {
        printf("FATAL ERROR: Unknown potential.\nABORTING SIMULATION\n\n");
        return NULL;
    }
    
    printMCP(mcs);
    
    
    // Number of pairs is calculated from total number of particles.
    //  At this time, arrays are set up to enable unlimited neighbor
    //  interactions.  In the future, these arrays could be smaller
    //  for systems with limited neighbors.
    mcs->numPairs   = ((mcs->N-1)*mcs->N)/2;
    
    // Set up arrays defining particles and pairs
    mcs->r          = (double *) malloc(mcs->N*sizeof(double));
    mcs->rij        = (double *) malloc(mcs->numPairs*sizeof(double));
    
    // Initialize arrays defining trial particles and pairs
    mcs->rTrial     = (double *) malloc(mcs->N*sizeof(double));
    mcs->rijTrial   = (double *) malloc(mcs->numPairs*sizeof(double));

    // Initialize arrays used to double-check the system energy
    rijTest         = (double *) malloc(mcs->numPairs*sizeof(double));
    
    // Initialize arrays used to convert between particle identifiers, ii and jj, -
    //  and pair identifier, ind.
    mcs->iii        = (unsigned long int *) malloc(mcs->numPairs*sizeof(double));
    mcs->jjj        = (unsigned long int *) malloc(mcs->numPairs*sizeof(double));
    mcs->indind     = (unsigned long int **) malloc((mcs->N-1)*sizeof(unsigned long int *));
    
    // Initialize the array ("indind") that maps the pair of particle indices
    //   ii and dj (= jj-ii-1) to a scalar pair index.
    // Also initialize the "inverse" arrays ("iii","jjj") that map the scalar pair index
    //   back to the pair of scalar particle indices.
    ind = 0;
    for (ii = 0; ii < mcs->N-1; ii++) {
            (*mcs).indind[ii] = (unsigned long int *) malloc((mcs->N-1-ii)*sizeof(unsigned long int));
            for (jj = ii + 1; jj < mcs->N; jj++) {
                    dj = jj - ii - 1;
                    mcs->indind[ii][dj] = ind;
                    mcs->iii[ind] = ii;
                    mcs->jjj[ind] = jj;
                    ind++;
            }
    }
    
    // Initialize accumulators to 0
    mcs->rhoSysA    = 0;
    mcs->rhoSysSA   = 0;
    mcs->lA         = 0;
    mcs->lSA        = 0;
    mcs->EA         = 0;
    mcs->ESA        = 0;
    mcs->lEA        = 0;
    mcs->VirA       = 0;
    mcs->VirSA      = 0;

    mcs->dAcc[0]    = 0;
    mcs->dAcc[1]    = 0;
    mcs->vAcc[0]    = 0;
    mcs->vAcc[1]    = 0;
    
    // Initialize block accumulator, last-step-state, and block mean arrays
    //  for density.
    mcs->rhol = (int *) malloc(mcs->rhonb*sizeof(int));
    mcs->rhoA = (int *) malloc(mcs->rhonb*sizeof(int));
    mcs->rhoM = (double *) malloc(mcs->rhonb*sizeof(double));
    for (ii = 0; ii < mcs->rhonb; ii++) {
            mcs->rhol[ii] = 0;
            mcs->rhoA[ii] = 0;
            mcs->rhoM[ii] = 0;
    }
        
    // Initialize block accumulator, last-step-state, and block mean arrays
    //  for g(r) or two-particle density.  And, open the corresponding
    //  output files.
    mcs->gl = (int **) malloc(mcs->gns*sizeof(int *));
    mcs->gA = (int **) malloc(mcs->gns*sizeof(int *));
    mcs->gM = (double **) malloc(mcs->gns*sizeof(double *));
    mcs->gf = (FILE **) malloc(mcs->gns*sizeof(FILE *));
    for (ii = 0; ii < mcs->gns; ii++) {
                mcs->gl[ii] = (int *) malloc(mcs->gnb*sizeof(int));
                mcs->gA[ii] = (int *) malloc(mcs->gnb*sizeof(int));
                mcs->gM[ii] = (double *) malloc(mcs->gnb*sizeof(double));
                for (jj = 0; jj < mcs->gnb; jj++) {
                        mcs->gl[ii][jj] = 0;
                        mcs->gA[ii][jj] = 0;
                        mcs->gM[ii][jj] = 0;
                } 
    }
    
    // Initialize values for a fresh (non-restart) run
    if ( mcs->isRestart == false ) {
        
        // Open output files
        mcs->cf   = fopen("config.dat.mcs","w");
        mcs->tf   = fopen("thermo.dat.mcs","w");
        mcs->rhof = fopen("rho.dat.mcs","w");
        
        for (ii = 0; ii < mcs->gns; ii++) {
                sprintf(gfstr,"g%lu.dat.mcs",ii);
                mcs->gf[ii] = fopen(gfstr,"w");
        }
        
        mcs->sn         = 0;
        mcs->slcp       = -1; 
        mcs->sltp       = -1;
        mcs->slrho      = -1;
        mcs->slg        = -1;
        
        // Initialize energy and virial terms to large values
        mcs->E          = 10E10;
        mcs->Vir        = 10E10;
        mcs->md         = 0;
        
        // Initialize box size based on number of particles.
        mcs->l          = mcs->N;

        // Initial values for particle positions are calculated by
        //  the distributing the particles evenly in the box.
        for (ii = 0; ii < mcs->N; ii++) {
            mcs->r[ii]        = -mcs->l/2 + (ii+0.5)*(mcs->l/mcs->N);
        }
        
        // Print headers to output files.
        fprintf(mcs->tf,"Step     Econf  Econf2   L     L2     LEconf  rho     rho2    Virial  Virial2\n");
        
    }
    
    // If this is a restart run, read info from data files to re-initialize
    else if ( mcs->isRestart == true ) {
        
        bool lineFound,stepFound,completeConfig=false;
        unsigned int pn,lineSize = 20000;
        unsigned long int lastrho,lastrhotmp,lastthermo,lastthermotmp,lastg[mcs->gns],lastgtmp;
        int nc,linesBelow=-2,npexp=0;
        long int cfpp = 0,tfpp = 0,rhofpp = 0,gfpp[mcs->gns]; // config file pointer position, etc.
        double x,y,z;
        char line1[lineSize],line2[lineSize];
        
        // The final time step in the configuration file drives the whole restart process
        //   because you cannot get back on the same track unless you know the 
        //   configuration at the point you want to re-start.
        //
        // *---------- RE-START CONFIGURATION FILE -----------*
        printf("Opening file %s\n","config.dat.mcs");
        mcs->cf = fopen("config.dat.mcs","r+");
        // Find last complete configuration in config.dat.mcs.
        // Read data from that configuration.
        
        // Setup the step number and particle positions from config.mcs.dat
        //   Start at the end of the file
        fseek(mcs->cf,-2,SEEK_END);
        
        // Continue looping until a COMPLETE configuration is found (all particles)
        while ( !completeConfig ) {
            stepFound = false;
            
            // Loop backwards over lines until 'Step' is found
            while ( !stepFound ) {
                linesBelow++;  // Lines counted (in reverse) from EOF or last incomplete 'Step'
                lineFound = false;
                
                // Store 2 lines at a time, makes it easier to get the number of particles
                //   expected, once 'Step' is found
                strncpy(line2,line1,lineSize);
                strncpy(line1,"",lineSize);

                nc = 0;  // Counter for the number of characters on a line
                
                // Loop backwards through characters until a new line is found
                while ( !lineFound ) {
                    nc++;
                    if ( ftell(mcs->cf) == 0) {
                        // Reached beginning of file
                        lineFound = true;
                        fgets(line1, lineSize, mcs->cf);
                        fseek(mcs->cf,0,SEEK_SET);
                        cfpp = ftell(mcs->cf);
                    }
                    else if ( fgetc(mcs->cf) == '\n') {
                        // Reached beginning of the line
                        lineFound = true;
                        fgets(line1, lineSize, mcs->cf);
                        fseek(mcs->cf,-nc-2,SEEK_CUR);
                        cfpp = ftell(mcs->cf) + 2;
                    }
                    else {
                        // Have not reached the beginning of the line yet
                        fseek(mcs->cf,-2,SEEK_CUR);
                    }
                }
                
                // If the LAST line (not this line) contained the substring 'Step'
                //   then read the step number and box length from that line
                //   Notice that the step number and box length will be overwritten
                //   later, if this step is found to be incomplete.
                if ( strstr(line2, "Step") ) {
                    stepFound = true;
                    sscanf(line2,"%*s %*s %lu %*s %*s %lg",&(mcs->sn),&(mcs->l));
                }
            }
            
            sscanf(line1,"%d",&npexp); // Read the expected number of particles from this line
            
            // If the expected number of particles is not the same as in INPUT, something
            //   is not right.  Issue an error, and exit.
            if ( npexp != mcs->N ) {
                printf("ERROR: Number of particles in config.mcs.dat (%d) is incompatible with "
                       "number of particles in INPUT (%lu).\n",npexp,mcs->N);
                exit(1);
            }
            
            // If the number of lines counted up from below matches the number of particles
            //   in the system, then this is a complete configuration.  Exit the loop.
            if ( linesBelow == mcs->N) {
                printf("Complete configuration found for step number %lu\n",mcs->sn);
                printf("Box length: %.8g\n",mcs->l);
                completeConfig = true;
            }
            // If this is not a complete configuration, move the cursor back to the end
            //   of the previous line, and start looping backwards through the prior
            //   configuration.
            else {
                linesBelow = -2;
            }
        }
        
        // Move the cursor back to the line stating the number of particles in the configuration.
        //   Then skip over that line, and the subsequent line (step no. and box length values)
        //   to put the cursor at the start of the line containing the position of particle 1.
        fseek(mcs->cf,cfpp,SEEK_SET);
        fgets(line1,lineSize,mcs->cf);
        fgets(line1,lineSize,mcs->cf);
        
        // Initialize particle positions from the identified last complete configuration.
        for (ii = 0; ii < mcs->N; ii++) {
            fgets(line1,lineSize,mcs->cf);
            sscanf(line1,"%u %lg %lg %lg",&pn,&x,&y,&z);
            mcs->r[ii]        = z;
        }
        
        // Overwrite any subsequent incomplete configuration data with whitespace, and
        //   move the pointer back to its previous position at the end of the last 
        //   complete config
        long int pfcc = ftell(mcs->cf); // position of (end of) final complete config
        for (ii = 0; ii < mcs->N; ii++ ) {
            fprintf(mcs->cf,"\n");
        }
        fseek(mcs->cf,pfcc,SEEK_SET);
        mcs->slcp = mcs->sn; 



        // *---------- RE-START DENSITY FILE -----------*
        printf("Opening file %s\n","rho.dat.mcs");
        mcs->rhof = fopen("rho.dat.mcs","r+");
        // Put the cursor at the last step equal or prior to the step found from the config file.
        //   Data points beyond that will be overwritten.
        while ( fgets(line1,lineSize,mcs->rhof) ) {
            sscanf(line1,"%lu",&lastrhotmp);
            if ( lastrhotmp <= mcs->sn) {
                lastrho = lastrhotmp;
                rhofpp = ftell(mcs->rhof);
            }
        }
        fseek(mcs->rhof,rhofpp,SEEK_SET);
        printf("Last rho time point <= last complete config time point: %lu\n"
               "  Overwriting density data beyond this point.\n",lastrho);
        mcs->slrho = mcs->sn;  // mcs->slrho is used when the density data is printed
                               //   to properly average over the block
        
        
        // *----------- RE-START g(r) FILES ------------*
        for (ii = 0; ii < mcs->gns; ii++) {
            sprintf(gfstr,"g%lu.dat.mcs",ii);
            printf("Opening file %s\n",gfstr);
            mcs->gf[ii] = fopen(gfstr,"r+");
            // printf("Cursor position: %ld\n",ftell(mcs->gf[ii]));
        }
        
        // For each g(r)-file, put the cursor at the last step equal or prior to the
        //   step found from the config file.
        unsigned long int lastgmin = mcs->sn;
        for (ii = 0; ii < (unsigned long int) mcs->gns; ii++) {
            lastgtmp = 0;
            while ( fgets(line1,lineSize,mcs->gf[ii]) ) {
                sscanf(line1,"%lu %*g ",&(lastgtmp));
                if ( lastgtmp <= mcs->sn ) {
                    lastg[ii] = lastgtmp;
                    gfpp[ii]  = ftell(mcs->gf[ii]);
                }
            }
            fseek(mcs->gf[ii],gfpp[ii],SEEK_SET);
            if (lastg[ii] < lastgmin) lastgmin = lastg[ii];
            printf("Last g[%lu] time point <= last complete config time point: %lu\n",ii,lastg[ii]);
        }
        mcs->slg = mcs->sn;  // mcs->slg is used when the density data is printed
                             //   to properly average over the block
        printf("Least compatible g time point: %lu\n"
               "  Overwriting g(r) data beyond this point.\n",lastgmin);
        
        
        // *---------- RE-START THERMO FILE ------------*
        printf("Opening file %s\n","thermo.dat.mcs");
        mcs->tf = fopen("thermo.dat.mcs","r+");
        while ( lastthermo <= mcs->sn ) {
            tfpp = ftell(mcs->tf);
            fgets(line1,lineSize,mcs->tf);
            sscanf(line1,"%lu",&lastthermo);
        }
        fseek(mcs->tf,tfpp,SEEK_SET);
        printf("Last thermo time point <= last complete config time point: %lu\n"
               "  Overwriting thermo data beyond this point.\n",lastthermo);
        mcs->sltp = mcs->sn;
        
        
        // Initialize energy and virial terms from (re-)starting configuration.
        mcs->E          = 10.0E10;
        mcs->Vir        = 10.0E10;
        mcs->md         = 0;
        
    }
    
    for (ind = 0; ind < mcs->numPairs; ind++) {
            mcs->rij[ind] = mcs->r[mcs->jjj[ind]] - mcs->r[mcs->iii[ind]];
    }
    
    // Perform a full calculation of the current density and g(r) or two-
    //  particle density histogram.
    fgrho(mcs);
    
    // Update the density and g(r) accumulators.
    ugrho(mcs);
    
    // Create the random number generator for the system.
    mcs->rangen = gsl_rng_alloc(gsl_rng_taus2);
    // Seed the random generator 
    gsl_rng_set(mcs->rangen,mcs->seed);
    
    // Flush streams to ensure everything prints before a crash
    fflush(mcs->cf);
    fflush(mcs->tf);
    fflush(stdout);
    
    // Return the pointer to the MCState struct that has just been setup
    return mcs;
}

void freeMCS(struct MCState *mcs) {
	int ii;
	
	free(mcs->eij);
	free(mcs->e12ij);
	free(mcs->e6ij);
	free(mcs->virij);
	free(mcs->vir12ij);
	free(mcs->vir6ij);
	free(mcs->eijTrial);
	free(mcs->e12ijTrial);
	free(mcs->e6ijTrial);
	free(mcs->virijTrial);
	free(mcs->vir12ijTrial);
	free(mcs->vir6ijTrial);
	free(mcs->r);
	free(mcs->rij);
	free(mcs->rTrial);
	free(mcs->rijTrial);
	free(mcs->indind[ii]);
        for (ii = 0; ii < mcs->N-1; ii++) {
                free(mcs->indind[ii]);
        }
        free(mcs->indind);
	free(mcs->iii);
	free(mcs->jjj);
	free(mcs->rhol);
	free(mcs->rhoA);
	free(mcs->rhoM);
        for (ii = 0; ii < mcs->gns; ii++) {
                free(mcs->gl[ii]);
                free(mcs->gA[ii]);
                free(mcs->gM[ii]);
                fclose(mcs->gf[ii]);
	}
	free(mcs->gl);
	free(mcs->gA);
	free(mcs->gM);
	free(mcs->gf);
	free(rijTest);
	free(eijTest);
	free(e12ijTest);
	free(e6ijTest);
	
	free(mcs);
}

// Print the current step
void printStep(struct MCState *mcs) {
    printf("Step: %lu...\n",mcs->sn);
    fflush(stdout);
}


// Do a full attempted trial displacement (no short-cuts to calculate the
//  energy using knowledge of the previous state).
int fad(struct MCState *mcs,unsigned long int *nm,double *d) {
	unsigned long int ii,jj,ind;
	double rij1,rij3,rij6,rij7,rij12,rij13;
	double phiij[6],params[2];
	params[0] = 1;

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
	
	mcs->md = (mcs->rn - 0.5)*2*mcs->maxStep;
	mcs->rTrial[mcs->nm] = mcs->r[mcs->nm] + mcs->md;
	if (fabs(mcs->rTrial[mcs->nm]) > (mcs->l/2.0)) {
		mcs->bndchk = 1;
		mcs->dAcc[1]++;
		mcs->ETrial = 10E10;
	}
	else {
		mcs->bndchk = 0;
	}
	fflush(stdout);
	}
	
	mcsETrial      = 0;
	mcsE12Trial    = 0;
	mcsE6Trial     = 0;
	mcsVirTrial    = 0;
	mcsVir12Trial  = 0;
	mcsVir6Trial   = 0;
	
	if (mcs->bndchk == 0) {
		params[1] = mcs->l;
		#pragma omp for private(ii,jj,ind,rij1,rij3,rij6,rij7,rij12,rij13) \
				reduction(+:mcsETrial,mcsE6Trial,mcsE12Trial, \
                                   mcsVirTrial,mcsVir6Trial,mcsVir12Trial)
		for (ind = 0; ind < mcs->numPairs; ind++) {
			ii            = mcs->iii[ind];
			jj            = mcs->jjj[ind];
			mcs->rijTrial[ind] = mcs->rTrial[jj] - mcs->rTrial[ii];
			if (mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
				mcs->eijTrial[ind]      = 0;
				mcs->e6ijTrial[ind]     = 0;
				mcs->e12ijTrial[ind]    = 0;
				mcs->virijTrial[ind]    = 0;
				mcs->vir6ijTrial[ind]   = 0;
				mcs->vir12ijTrial[ind]  = 0;	
			}
			else {
				rij1          = mcs->rijTrial[ind];
				
				mcs->phi(&rij1,mcs->potCutOff,params,phiij);
				mcs->eijTrial[ind]     = phiij[0];
				mcs->e12ijTrial[ind]   = phiij[2];
                                mcs->e6ijTrial[ind]    = phiij[4];
                                mcs->virijTrial[ind]   = phiij[1];
				mcs->vir12ijTrial[ind] = phiij[3];
				mcs->vir6ijTrial[ind]  = phiij[5];
				
				mcsETrial      += mcs->eijTrial[ind];
				mcsE12Trial    += mcs->e12ijTrial[ind];
				mcsE6Trial     += mcs->e6ijTrial[ind];
				mcsVir12Trial  += mcs->vir12ijTrial[ind];
				mcsVirTrial    += mcs->virijTrial[ind];
				mcsVir6Trial   += mcs->vir6ijTrial[ind];
			}
		}
		
		mcs->ETrial      = mcsETrial;
		mcs->E12Trial    = mcsE12Trial;
		mcs->E6Trial     = mcsE6Trial;
		mcs->VirTrial    = mcsVirTrial;
		mcs->Vir12Trial  = mcsVir12Trial;
		mcs->Vir6Trial   = mcsVir6Trial;
		
		
		#pragma omp single
		{

		if (mcs->ETrial > mcs->E ) {
			ran = gsl_rng_uniform(mcs->rangen);
			mcs->bf = exp((mcs->E-mcs->ETrial)/mcs->T);
		}

		if (mcs->ETrial < mcs->E || mcs->bf > ran) {
			mcs->dAcc[0]++;
			mcs->E     = mcs->ETrial;
			mcs->E12   = mcs->E12Trial;
			mcs->E6    = mcs->E6Trial;
			mcs->Vir   = mcs->VirTrial;
			mcs->Vir12 = mcs->Vir12Trial;
			mcs->Vir6  = mcs->Vir6Trial;
			
			mcs->r[mcs->nm] = mcs->rTrial[mcs->nm];
			
			for (ind = 0; ind < mcs->numPairs; ind++) {
				mcs->rij[ind]     = mcs->rijTrial[ind];
				mcs->eij[ind]     = mcs->eijTrial[ind];
				mcs->e12ij[ind]   = mcs->e12ijTrial[ind];
				mcs->e6ij[ind]    = mcs->e6ijTrial[ind];
				mcs->virij[ind]   = mcs->virijTrial[ind];
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
				mcs->vir6ij[ind]  = mcs->vir6ijTrial[ind];
			}
		}
		else {
			mcs->dAcc[1]++;
		}
		}
	}
        //printf("After fad, E = %.5G\n",mcs->E);
	fflush(stdout);
	return 0;
	#pragma omp barrier
}


// Print coordinates of all particles in the system
int printCoords(struct MCState *mcs) {
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


// Print the density histogram
int printRho(struct MCState *mcs) {

        unsigned long int rb,ns;
        ns = (mcs->sn - mcs->slrho);
        fprintf(mcs->rhof,"%lu",mcs->sn);
        
        for (rb = 0; rb < mcs->rhonb; rb++) {
          mcs->rhoM[rb] = (double) mcs->rhoA[rb]/ns/mcs->rbw;
          fprintf(mcs->rhof," %.8G",mcs->rhoM[rb]);
          mcs->rhoA[rb] = 0;
        }
        
        fprintf(mcs->rhof,"\n");
        fflush(mcs->rhof);
        mcs->slrho = mcs->sn;

        return 0;
}


// Print the g(x) histograms
int printG(struct MCState *mcs) {
  
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
  
  
// Do a full calculation of density and g(x), not using any knowledge of the
//  distributions at the previous step.
int fgrho(struct MCState *mcs) {

	unsigned long int ii,jj,ind;
	long int rb,gs1,gs2,gb;
	
	#pragma omp for
	for (rb = 0; rb < mcs->rhonb; rb++) {
		mcs->rhol[rb] = 0;
	}

	#pragma omp for
	for(ii = 0; ii < mcs->N; ii++) {
		rb = (long int) floor(mcs->r[ii]/mcs->rbw + mcs->rhonb/2.0);
		if (rb >= 0 && rb < mcs->rhonb) {
			#pragma omp atomic
			mcs->rhol[rb]++;
		}
	}
	
	#pragma omp barrier
	
	#pragma omp for
	for (gs1 = 0; gs1 < mcs->gns; gs1++) {
		for (gb = 0; gb < mcs->gnb; gb++) {
			mcs->gl[gs1][gb] = 0;
		}
	}
	
	#pragma omp for private(ind,ii,jj,gs1,gs2,gb)
	for (ind = 0; ind < mcs->numPairs; ind++) {
		ii = mcs->iii[ind];
		jj = mcs->jjj[ind];
		gs1 = (long int) floor(mcs->r[ii]/mcs->gsw + mcs->gns/2.0);
		gs2 = (long int) floor(mcs->r[jj]/mcs->gsw + mcs->gns/2.0);
		gb = (long int) floor(mcs->rij[ind]/mcs->gbw);
	        
		if(gb < mcs->gnb) {
			if (gs1 >= 0 && gs1 < (long int) mcs->gns) {
				#pragma omp atomic
				mcs->gl[gs1][gb]++;
			}
			else {
				
			}
			
			if (gs2 >= 0 && gs2 < (long int) mcs->gns) {
				#pragma omp atomic
				mcs->gl[gs2][gb]++;
			}
			else {
			}
		
		}
		else {
		}
	}

  return 0;
}


// Accumulate density and g(x)
int ugrho(struct MCState *mcs) {

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


// Quickly attempt a displacement trial, using knowledge of the state at the
//  previous step. Update the state according to the outcome of the trial.
int qad(struct MCState *mcs,unsigned long int *nm, double *d) {
        int ii = mcs->qad(mcs,nm,d);
        return ii;
}


int qad2(struct MCState *mcs,unsigned long int *nm, double *d) {
        
	unsigned long int ii,jj,dj,ind;
	double rij1,rij3,rij6,rij7,rij12,rij13;
        double phiij[6],params[2];
	params[0] = 1;
	params[1] = mcs->l;

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
	mcs->md = (mcs->rn - 0.5)*2*mcs->maxStep;
	mcs->rTrial[mcs->nm] = mcs->r[mcs->nm] + mcs->md;
        //printf("nm: %lu  rn: %.5G  md: %.5G  rT: %.5G  \n",mcs->nm,mcs->rn,mcs->md,mcs->rTrial[mcs->nm]);
        //fflush(stdout);
	}
	
	if (fabs(mcs->rTrial[mcs->nm]) > mcs->l/2.0) {
		#pragma omp single
		{
		mcs->dAcc[1]++;
		}
	}
	else {
		//#pragma omp single
		{
		mcsdE = 0;
		mcsdE12 = 0;
		mcsdE6 = 0;
		mcsdVir = 0;
		mcsdVir12 = 0;
		mcsdVir6 = 0;
                //printf("mcs->dVir,mcsdVir = %5.5G, %5.5G\n",mcs->dVir,mcsdVir);
		}
		//#pragma omp barrier
                
		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) \
                            reduction(+:mcsdE,mcsdE6,mcsdE12,mcsdVir,mcsdVir6,mcsdVir12)
		for (ii = 0; ii < mcs->nm; ii++) {
			dj = mcs->nm - ii - 1;
			ind = mcs->indind[ii][dj];
			mcs->rijTrial[ind] = mcs->rij[ind] + mcs->md;
			if ( mcs->nbn >= 0 && (unsigned int) dj >= mcs->nbn ) {
				mcs->eijTrial[ind] = 0;
				mcs->e12ijTrial[ind] = 0;
                                mcs->e6ijTrial[ind] = 0;
				mcs->virijTrial[ind] = 0;
				mcs->vir12ijTrial[ind] = 0;
				mcs->vir6ijTrial[ind] = 0;
			}
			else {
				rij1 = mcs->rijTrial[ind];
				
				mcs->phi(&rij1,mcs->potCutOff,params,phiij);
				mcs->eijTrial[ind]     = phiij[0];
				mcs->e12ijTrial[ind]   = phiij[2];
                                mcs->e6ijTrial[ind]    = phiij[4];
                                mcs->virijTrial[ind]   = phiij[1];
				mcs->vir12ijTrial[ind] = phiij[3];
				mcs->vir6ijTrial[ind]  = phiij[5];
			        //printf("NEIGHBOR!!!! ii: %lu  dj: %lu  rij1: %.5G  eijTrial: %.5G  \n",ii,dj,rij1,mcs->eijTrial[ind]);
                                //fflush(stdout);
			}
			mcsdE     = mcsdE     - mcs->eij[ind]     + mcs->eijTrial[ind];
			mcsdE12   = mcsdE12   - mcs->e12ij[ind]   + mcs->e12ijTrial[ind];
			mcsdE6    = mcsdE6    - mcs->e6ij[ind]    + mcs->e6ijTrial[ind];
			/*if ( isnan(mcs->virijTrial[ind]) || isnan(mcs->virij[ind]) || isnan(mcsdVir) ) {
			//if ( mcs->sn == 6166) {
				printf("For ind = %lu,  r[%lu] = %5.5G, r[%lu] = %5.5G, rij,(trial) = %5.5G,%5.5G, virij,(trial) = %5.5G,%5.5G, mcsdVir = %5.5G\n",ind,ii,mcs->r[ii],mcs->nm,mcs->r[mcs->nm],mcs->rij[ind],mcs->rijTrial[ind],mcs->virij[ind],mcs->virijTrial[ind],mcsdVir);
			}*/
			/*if ( isnan(mcs->virij[ind])) {
			//if ( mcs->sn == 6166) {
				printf("For ind = %lu,  r[%lu] = %5.5G,  mcs->virij = %5.5G!\n",ind,ii,ii+dj+1,mcs->virij[ind]);
			}
			if ( isnan(mcsdVir)) {
			//if ( mcs->sn == 6166) {
				printf("For ind = %lu,  r[%lu][%lu] = %5.5G,  mcsdVir (before adding ii,jj contributions) = %5.5G!\n",ii,ii+dj+1,ind,mcsdVir);
			}*/
			mcsdVir   = mcsdVir   - mcs->virij[ind]   + mcs->virijTrial[ind];
			/*if ( isnan(mcsdVir) ) {
				printf("For ii,jj,ind = %lu,%lu,%lu,  mcsdVir (after adding ii,jj contributions) = %5.5G!\n",ii,ii+dj+1,ind,mcsdVir);
			}*/
			mcsdVir12 = mcsdVir12 - mcs->vir12ij[ind] + mcs->vir12ijTrial[ind];
			mcsdVir6  = mcsdVir6  - mcs->vir6ij[ind]  + mcs->vir6ijTrial[ind];
		}
                //printf("mcsdVir = %5.5G\n",mcsdVir);
                //fflush(stdout);
                #pragma omp barrier
		#pragma omp single
		{
		//if ( isnan(mcsdVir) ) {
                //    printf("dVir is NAN!!\n");
                //}
                mcs->dE = mcsdE;
                mcs->dE12 = mcsdE12;
		mcs->dE6 = mcsdE6;
		mcs->dVir = mcsdVir;
		mcs->dVir12 = mcsdVir12;
		mcs->dVir6 = mcsdVir6;
		//printf("mcs->dVir,mcsdVir: %5.5G, %5.5G\n",mcs->dVir,mcsdVir);
                //fflush(stdout);
		}

                {
                mcsdE = 0;
		mcsdE12 = 0;
		mcsdE6 = 0;
		mcsdVir = 0;
		mcsdVir12 = 0;
		mcsdVir6 = 0;
                }
		
                ii = 7;
                //printf("mcs->nm: %lu  ii: %lu\n",mcs->nm,ii);
                //fflush(stdout);
		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) \
                            reduction(+:mcsdE,mcsdE6,mcsdE12,mcsdVir,mcsdVir6,mcsdVir12)
		for (ii = (mcs->nm + 1); ii < mcs->N; ii++) {
                        dj = ii - mcs->nm - 1;
			ind = mcs->indind[mcs->nm][dj];
			mcs->rijTrial[ind] = mcs->rij[ind] - mcs->md;
			if ( mcs->nbn >= 0 && (unsigned int) dj >= mcs->nbn ) {
				mcs->eijTrial[ind] = 0;
                                mcs->e6ijTrial[ind] = 0;
				mcs->e12ijTrial[ind] = 0;
				mcs->virijTrial[ind] = 0;
				mcs->vir6ijTrial[ind] = 0;
				mcs->vir12ijTrial[ind] = 0;
			}
			else {
				rij1 = mcs->rijTrial[ind];
				
				mcs->phi(&rij1,mcs->potCutOff,params,phiij);
				mcs->eijTrial[ind]     = phiij[0];
				mcs->e12ijTrial[ind]   = phiij[2];
                                mcs->e6ijTrial[ind]    = phiij[4];
                                mcs->virijTrial[ind]   = phiij[1];
				mcs->vir12ijTrial[ind] = phiij[3];
				mcs->vir6ijTrial[ind]  = phiij[5];
			        //printf("NEIGHBOR!!!! ii: %lu  dj: %lu  rij1: %.5G  eijTrial: %.5G  \n",ii,dj,rij1,mcs->eijTrial[ind]);
                                //fflush(stdout);
                        }
			mcsdE     = mcsdE     - mcs->eij[ind]     + mcs->eijTrial[ind];
                        mcsdE12   = mcsdE12   - mcs->e12ij[ind]   + mcs->e12ijTrial[ind];
			mcsdE6    = mcsdE6    - mcs->e6ij[ind]    + mcs->e6ijTrial[ind];
			mcsdVir   = mcsdVir   - mcs->virij[ind]   + mcs->virijTrial[ind];
			mcsdVir12 = mcsdVir12 - mcs->vir12ij[ind] + mcs->vir12ijTrial[ind];
			mcsdVir6  = mcsdVir6  - mcs->vir6ij[ind]  + mcs->vir6ijTrial[ind];
		}
                //printf("mcsdVir = %5.5G\n",mcsdVir);
                //fflush(stdout);
		#pragma omp barrier  // THIS BARRIER IS CRITICAL TO CALCULATE dE ACCURATELY (I DO NOT UNDERSTAND WHY THAT SHOULD BE THE CASE)
		#pragma omp single
		{
		mcs->dE = mcs->dE + mcsdE;
                mcs->dE12 = mcs->dE12 + mcsdE12;
		mcs->dE6 = mcs->dE6 + mcsdE6;
		mcs->dVir = mcs->dVir + mcsdVir;
		mcs->dVir12 = mcs->dVir12 + mcsdVir12;
		mcs->dVir6 = mcs->dVir6 + mcsdVir6;
		//printf("mcs->dVir,mcsdVir: %5.5G, %5.5G\n",mcs->dVir,mcsdVir);
                //fflush(stdout);
                }
		
		if (mcs->dE > 0 ) {
		#pragma omp single
		{
			ran = gsl_rng_uniform(mcs->rangen);
			mcs->bf = exp(-mcs->dE/mcs->T);
		        //printf("ran: %.5G  bf: %.5G  ",ran,mcs->bf);
                        //fflush(stdout);
                }
		}
		
		if (mcs->dE <= 0 || mcs->bf > ran) {
			#pragma omp single
			{
			//printf("At mid-qad (%lu %.5G->(%.5G)->%.5G), E = %.7G, dE = %.7G\nE should = %.5G\n",
                        //         mcs->nm,mcs->r[mcs->nm],mcs->md,mcs->rTrial[mcs->nm],mcs->E,
                        //         mcsdE,mcs->E + mcsdE);
			mcs->dAcc[0]++;
			mcs->E      += mcs->dE;
			mcs->E12    += mcs->dE12;
			mcs->E6     += mcs->dE6;
			mcs->Vir    += mcs->dVir;
			mcs->Vir12  += mcs->dVir12;
			mcs->Vir6   += mcs->dVir6;
			mcs->r[mcs->nm] = mcs->rTrial[mcs->nm];
			//printf("Accepted  ");
                        //fflush(stdout);
			}
			
			#pragma omp for private(jj,dj,ind) nowait
			for (jj = 0; jj < mcs->nm; jj++) {
				dj = mcs->nm - jj - 1;
				ind = mcs->indind[jj][dj];
				mcs->rij[ind] = mcs->rijTrial[ind];
				mcs->eij[ind] = mcs->eijTrial[ind];
				mcs->e12ij[ind] = mcs->e12ijTrial[ind];
				mcs->e6ij[ind] = mcs->e6ijTrial[ind];
				mcs->virij[ind] = mcs->virijTrial[ind];
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
				mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
			}
			#pragma omp for private(jj,dj,ind)
			for (jj = mcs->nm + 1; jj < mcs->N; jj++) {
				dj = jj - mcs->nm - 1;
				ind = mcs->indind[mcs->nm][dj];
				mcs->rij[ind] = mcs->rijTrial[ind];
				mcs->eij[ind] = mcs->eijTrial[ind];
				mcs->e12ij[ind] = mcs->e12ijTrial[ind];
				mcs->e6ij[ind] = mcs->e6ijTrial[ind];
				mcs->virij[ind] = mcs->virijTrial[ind];
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
				mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
			}
			
                        // Short-cut calculation of density and g(x)
			qagrho(mcs);
			/*
			#pragma omp single
			{
			printf("r[28]: %.5G  r[32]: %.5G  rij[28][32]: %.5G\n",mcs->r[28],mcs->r[32],mcs->rij[mcs->indind[28][3]]);
			printf("r[32]: %.5G  r[34]: %.5G  rij[32][34]: %.5G\n",mcs->r[32],mcs->r[34],mcs->rij[mcs->indind[32][1]]);
			printf("r[32]: %.5G  r[35]: %.5G  rij[32][35]: %.5G\n",mcs->r[32],mcs->r[35],mcs->rij[mcs->indind[32][2]]);
			fflush(stdout);
			}
			*/
		}
		else {
		#pragma omp single
		{
                        //printf("Rejected  ");
                        //fflush(stdout);
			mcs->dAcc[1]++;
		}
		}
	}
	
        // Accumulate density and g(x)
	ugrho(mcs);
	
        #pragma omp single
        {
          //printf("mcs->E: %.5G\n",mcs->E);
          //fflush(stdout);
        }
        fflush(stdout);
	#pragma omp barrier

	return 0;
}


/*
int qadLJ(struct MCState *mcs,unsigned long int *nm, double *d) {
        
	unsigned long int ii,jj,dj,ind;
	double rij1,rij3,rij6,rij7,rij12,rij13;
        
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

	mcs->md = (mcs->rn - 0.5)*2*mcs->maxStep;
	mcs->rTrial[mcs->nm] = mcs->r[mcs->nm] + mcs->md;
	}
	
	if (fabs(mcs->rTrial[mcs->nm]) > mcs->l/2.0) {
		#pragma omp single
		{
		mcs->dAcc[1]++;
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
		
		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) \
                    reduction(+:mcsdE,mcsdE6,mcsdE12,mcsdVir,mcsdVir6,mcsdVir12) nowait
		for (ii = 0; ii < mcs->nm; ii++) {
			dj = mcs->nm - ii - 1;
			ind = mcs->indind[ii][dj];
			mcs->rijTrial[ind] = mcs->rij[ind] + mcs->md;
			if ( mcs->nbn >= 0 && dj >= mcs->nbn ) {
				mcs->e6ijTrial[ind] = 0;
				mcs->e12ijTrial[ind] = 0;
				mcs->vir6ijTrial[ind] = 0;
				mcs->vir12ijTrial[ind] = 0;
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
			}
			mcsdE6  = mcsdE6  - mcs->e6ij[ind]  + mcs->e6ijTrial[ind];
			mcsdE12 = mcsdE12 - mcs->e12ij[ind] + mcs->e12ijTrial[ind];
			mcsdVir6 = mcsdVir6 - mcs->vir6ij[ind] + mcs->vir6ijTrial[ind];
			mcsdVir12 = mcsdVir12 - mcs->vir12ij[ind] + mcs->vir12ijTrial[ind];
		}

		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) \
                      reduction(+:mcsdE,mcsdE6,mcsdE12,mcsdVir,mcsdVir6,mcsdVir12)
		for (ii = mcs->nm + 1; ii < mcs->N; ii++) {
			dj = ii - mcs->nm - 1;
			ind = mcs->indind[mcs->nm][dj];
			mcs->rijTrial[ind] = mcs->rij[ind] - mcs->md;
			if ( mcs->nbn >= 0 && dj >= mcs->nbn ) {
				mcs->e6ijTrial[ind]    = 0;
				mcs->e12ijTrial[ind]   = 0;
				mcs->vir6ijTrial[ind]  = 0;
				mcs->vir12ijTrial[ind] = 0;
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
			}
			mcsdE6  = mcsdE6  - mcs->e6ij[ind]  + mcs->e6ijTrial[ind];
			mcsdE12 = mcsdE12 - mcs->e12ij[ind] + mcs->e12ijTrial[ind];
			mcsdVir6 = mcsdVir6 - mcs->vir6ij[ind] + mcs->vir6ijTrial[ind];
			mcsdVir12 = mcsdVir12 - mcs->vir12ij[ind] + mcs->vir12ijTrial[ind];
		}
		
		#pragma omp barrier  // THIS BARRIER IS CRITICAL TO CALCULATE dE ACCURATELY (I DO NOT UNDERSTAND WHY THAT SHOULD BE THE CASE)
		
		#pragma omp single
		{
		mcs->dE12 = mcsdE12;
		mcs->dE6 = mcsdE6;
		mcs->dE = mcs->dE12 - mcs->dE6;
		}
		
		if (mcs->dE > 0 ) {
		#pragma omp single
		{
			ran = gsl_rng_uniform(mcs->rangen);
			mcs->bf = exp(-mcs->dE/mcs->T);
		}
		}
		
		if (mcs->dE <= 0 || mcs->bf > ran) {
			#pragma omp single
			{
			mcs->dAcc[0]++;
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
				dj = mcs->nm - jj - 1;
				ind = mcs->indind[jj][dj];
				mcs->rij[ind] = mcs->rijTrial[ind];
				mcs->e6ij[ind] = mcs->e6ijTrial[ind];
				mcs->e12ij[ind] = mcs->e12ijTrial[ind];
				mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
			}
			#pragma omp for private(jj,dj,ind)
			for (jj = mcs->nm + 1; jj < mcs->N; jj++) {
				dj = jj - mcs->nm - 1;
				ind = mcs->indind[mcs->nm][dj];
				mcs->rij[ind] = mcs->rijTrial[ind];
				mcs->e6ij[ind] = mcs->e6ijTrial[ind];
				mcs->e12ij[ind] = mcs->e12ijTrial[ind];
				mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
			}
			
                        // Short-cut calculation of density and g(x)
			qagrho(mcs);
			
		}
		else {
		#pragma omp single
		{
			mcs->dAcc[1]++;
		}
		}
	}
	
        // Accumulate density and g(x)
	ugrho(mcs);
	
        fflush(stdout);
	#pragma omp barrier

	return 0;
}
*/

// Quickly attempt a volume change trial, using knowledge of the state at the
//  previous step. Update the state according to the outcome of the trial.
int qav(struct MCState *mcs) {
        int ii = mcs->qav(mcs);
        return ii;
}

	
int qavLJ(struct MCState *mcs) {
	unsigned long int ii,jj,ind;
	double dQ;
	double phiij[6],params[2];
	params[0] = 1;

	//printf("Into qavLJ\n");
	#pragma omp single
	{
	mcs->dl = (mcs->rn-0.5)*2*mcs->maxdl;
	lRat1 = (mcs->l+mcs->dl)/mcs->l;
	lRat3 = lRat1*lRat1*lRat1;
	lRat6 = 1/(lRat3*lRat3);
	lRat12 = lRat6*lRat6;
	mcs->E12Trial = lRat12*mcs->E12;
	mcs->E6Trial = lRat6*mcs->E6;
	mcs->dE = mcs->E12Trial - mcs->E6Trial - mcs->E;
	mcs->bf   = exp(-(mcs->dE + mcs->P*mcs->dl)/mcs->T + mcs->N*log(lRat1));

	if (mcs->bf < 1.0 ) {
		ran = gsl_rng_uniform(mcs->rangen);
	}
	}
	
	if (mcs->bf >= 1.0 || mcs->bf > ran) {
		#pragma omp single
		{
		mcs->vAcc[0]++;
		mcs->E   = mcs->E + mcs->dE;
		mcs->E12 = mcs->E12Trial;
		mcs->E6  = mcs->E6Trial;
		
		mcs->l = mcs->l + mcs->dl;
		
		lRat7 = lRat6/lRat1;
		lRat13 = lRat12/lRat1;
		mcs->Vir6 = lRat7*mcs->Vir6;
		mcs->Vir12 = lRat13*mcs->Vir12;
		mcs->Vir    = mcs->N*mcs->T/mcs->l + mcs->Vir12 - mcs->Vir6;
		}
		params[1] = (double) mcs->l;
		
		#pragma omp single private(ii) nowait
		for (ii = 0; ii < mcs->N; ii++) {
			mcs->r[ii] = lRat1*mcs->r[ii];
		}
		
		#pragma omp for private(ind)
		for (ind = 0; ind < mcs->numPairs; ind++) { 
			ii = mcs->iii[ind];
			jj = mcs->jjj[ind];
			mcs->rij[ind]     = lRat1*mcs->rij[ind];
      			mcs->phi(&(mcs->rij[ind]),mcs->potCutOff,&params,phiij);
	      		mcs->eij[ind] = phiij[0];
      			mcs->e12ij[ind] = phiij[2];
      			mcs->e6ij[ind] = phiij[4];
			mcs->virij[ind]  = phiij[1];
			/*if ( isnan( mcs->virij[ind] ) ) {
				printf("14 mcs->dl = %5.5G, mcs->l = %5.5G, params[1] = %5.5G\n",mcs->dl,mcs->l,params[1]);
				printf("PARAMS IN JMMMCSTATE: params,&(params[0]),&(params[1]) = %p,%p,%p  *params = %2.2G,%5.5G\n",params,&(params[0]),&(params[1]),params[0],params[1]);
				printf("r[%lu],r[%lu] = %5.5G,%5.5G  rij = %5.5G   virij = %5.5G\n",ii,jj,mcs->r[ii],mcs->r[jj],mcs->rij[ind],mcs->virij[ind]);
			}*/
			mcs->vir12ij[ind] = phiij[3];
			mcs->vir6ij[ind]  = phiij[5];
		}
		
		fgrho(mcs);
	}
	
	else {
		#pragma omp single
		{
		mcs->vAcc[1]++;
		}
	}
	ugrho(mcs);
	fflush(stdout);

	return 0;
}


// Increment the state's current step
unsigned long int incrementStep(struct MCState *mcs) {
  int flag = 0;
  const int stepPrintInterval = 10000;
  
  if (mcs->sn == mcs->numSteps) {
    return 0;
  }
  else {
    #pragma omp barrier
    #pragma omp single
    {
    mcs->sn++;
    
    if (mcs->sn % stepPrintInterval == 0) {
        printStep(mcs);
    }
    }
    return mcs->sn;
  }
  #pragma omp barrier
}


// Attempt a Monte Carlo step on the state.
int Step(struct MCState *mcs) {
  #pragma omp barrier
  #pragma omp single
  {
  mcs->nm = gsl_rng_uniform_int(mcs->rangen,mcs->N+1);
  mcs->rn = gsl_rng_uniform(mcs->rangen);
  }
  #pragma omp barrier

  // nm is the number of the particle selected to try displacing. If its
  //  number exceeds that of the highest-numbered particle, that means
  //  do a volume change trial.
  if (mcs->nm < mcs->N) {
    // Do a trial displacement on particle nm
   //printf("Displacing!  ************ *   * *   *   *   *   *\n");
    qad(mcs,&(mcs->nm),&(mcs->rn));
  }
  else {
    // Do a volume change trial
    // printf("qav-ing!  ************ *   * *   *   *   *   *\n");
    qav(mcs);
  }
  
  // At appropriate steps, double-check that the energy has not been
  //  corrupted.  (Useful when qad and qav are being used because they
  //  calculate the energy based on the change in energy from the 
  //  previous step.  Therefore, if a step's energy is calculated
  //  incorrectly, all subsequent steps will be calculated incorrectly
  //  until corrected.)
  if (isECheck(mcs) ) {
    ECheck(mcs);
  }
  #pragma omp single
  {
  updateThermo(mcs);
  }
  fflush(stdout);
  
  return 0;

}


// Check if this is a step at which coordinates should be printed
int isCoordPrint(struct MCState *mcs) {
  if ( mcs->sn % mcs->cpi == 0 ) {
    return 1;
  }
  else {
    return 0;
  }
}

// Check if this is a step at which thermodynamic parameters should be printed
int isThermoPrint(struct MCState *mcs ) {
  if ( mcs->sn % mcs->tpi == 0 ) {
    return 1;
  }
  else {
    return 0;
  }
}


// Check if this is a step at which density should be printed
int isRhoPrint(struct MCState *mcs) {
  if ( mcs->sn % mcs->rhopi == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}

// Check if this is a step at which g(x) should be printed
int isGPrint(struct MCState *mcs) {
  if ( mcs->sn % mcs->gpi == 0 ) {
    return 1;
  }
  else {
    return 0;
  }
}


// Check if this is a step at which energy should be checked
int isECheck(struct MCState *mcs) {
  if ( mcs->sn % mcs->eci == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}

// Check if this is a step at which maximum trial displacement should be
//  updated.
int isMaxDisAdjust(struct MCState *mcs) {
  
  if ( mcs->sn % mcs->mdai == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}

// Check if this is a step at which maximum trial volume change should be
//  updated.
int isMaxDVAdjust(struct MCState *mcs) {
  
  if (mcs->sn % mcs->mvai == 0 ) {
    return 1;
  }
  else {
    return 0;
  }

}


// Print thermodynamic properties, averaged over thermo block
int printThermo(struct MCState *mcs) {
	double rhoSysM,rhoSysSM,lM,lSM,EM,ESM,lEM,VirM,VirSM;
	unsigned long int ss;

	ss = mcs->sn - mcs->sltp;

	rhoSysM        = mcs->rhoSysA/ss;
	rhoSysSM       = mcs->rhoSysSA/ss;
	lM             = mcs->lA/ss;
	lSM            = mcs->lSA/ss;
	EM             = mcs->EA/ss;
	ESM            = mcs->ESA/ss;
        lEM            = mcs->lEA/ss;
	VirM           = mcs->VirA/ss;
	VirSM          = mcs->VirSA/ss;
	
        fprintf(mcs->tf,"%lu\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\t%.8G\n",mcs->sn,EM,ESM,lM,lSM,lEM,rhoSysM,rhoSysSM,VirM,VirSM);
	fflush(mcs->tf);
        printf("%lu  %.8G  %.8G  %.8G\n",mcs->sn,mcs->E,mcs->l,mcs->Vir);
        fflush(stdout);
	mcs->rhoSysA   = 0;
	mcs->rhoSysSA  = 0;
	mcs->lA        = 0;
	mcs->lSA       = 0;
	mcs->EA        = 0;
	mcs->ESA       = 0;
        mcs->lEA       = 0;
	mcs->VirA      = 0;
	mcs->VirSA     = 0;
	mcs->sltp      = mcs->sn;

	return 0;
}


// Accumulate thermodynamic property values
int updateThermo(struct MCState *mcs) {
    
    double rhotmp;
   
    rhotmp        = mcs->N/mcs->l;
    mcs->rhoSysA  = mcs->rhoSysA + rhotmp;
    mcs->rhoSysSA = mcs->rhoSysSA + rhotmp*rhotmp; 
    mcs->lA       = mcs->lA + mcs->l;
    mcs->lSA      = mcs->lSA + mcs->l*mcs->l; 
    mcs->EA       = mcs->EA + mcs->E; 
    mcs->ESA      = mcs->ESA + mcs->E*mcs->E; 
    mcs->lEA      = mcs->lEA + mcs->l*mcs->E; 
    mcs->VirA     = mcs->VirA + mcs->Vir; 
    mcs->VirSA    = mcs->VirSA + mcs->Vir*mcs->Vir;

    return 0;

}


// Double-check the system energy
int ECheck(struct MCState *mcs) {
  unsigned long int iq,ind,ii,jj;
  double rij1,rij3,rij6,rij12;
  double phiij[6],params[2];
  params[0] = 0;
  
  #pragma omp barrier 
  mcsETest = 0;
  #pragma omp for private(iq,ii,jj,rij1,rij3,rij6,rij12,phiij) reduction(+:mcsETest)
  for (iq = 0; iq < mcs->numPairs; iq++) {
    ii = mcs->iii[iq];
    jj = mcs->jjj[iq];
    
    rijTest[iq] = mcs->r[jj] - mcs->r[ii];
    if (mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
      eijTest[iq] = 0;
      e12ijTest[iq] = 0;
      e6ijTest[iq]  = 0;
    }
    else {
      rij1 = rijTest[iq];
      
      mcs->phi(&rij1,mcs->potCutOff,params,phiij);
      eijTest[iq] = phiij[0];
      e12ijTest[iq] = phiij[2];
      e6ijTest[iq] = phiij[4];
      mcsETest += eijTest[iq];
    }
  }
  
  //printf("Final mcsETest: %.5G\n",mcsETest);
  #pragma omp single
  {
  if (fabs(mcsETest - mcs->E) > 0.0001) {
    resultFlag = 1;
    fprintf(stdout,"\n### !!! $$$ Energy discrepancy!!! ######### !!!!!!!!!!! $$$$$"
         "$$$$$ !!!!!!!!\nWhile attempting to move particle: %lu by %.5G\nE = %.8G, ETest = %.8G, Resetting results "
         "to that of full calculation (E = %.8G)\n\n",mcs->nm,mcs->md,mcs->E,mcsETest,mcsETest);

    params[0] = 1;
    params[1] = mcs->l;
    mcs->E = mcsETest;
    mcsETrial     = 0;
    mcsE12Trial   = 0;
    mcsE6Trial    = 0;
    mcsVirTrial   = 0;
    mcsVir12Trial = 0;
    mcsVir6Trial  = 0;
  }
  }

  if (resultFlag == 1) {
  #pragma omp single
  {
  printf("Looking for particles contributing to the discrepancy.\n");
  }

  #pragma omp for private(ii,jj,ind,rij1,phiij) \
      reduction (+:mcsETrial,mcsE12Trial,mcsE6Trial,mcsVirTrial,mcsVir12Trial,mcsVir6Trial)
  for (ind = 0; ind < mcs->numPairs; ind++) {
    ii            = mcs->iii[ind];
    jj            = mcs->jjj[ind];

    if ( mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
      mcs->eij[ind]        = 0;
      mcs->e12ij[ind]      = 0;
      mcs->e6ij[ind]       = 0;
      mcs->virij[ind]      = 0;
      mcs->vir12ij[ind]    = 0;
      mcs->vir6ij[ind]     = 0;
    }
    else {
      mcs->rij[ind]          = mcs->r[jj] - mcs->r[ii];
			
      mcs->phi(&(mcs->rij[ind]),mcs->potCutOff,params,phiij);
      mcs->eij[ind]      = phiij[0];
      mcs->e12ij[ind]    = phiij[2];
      mcs->e6ij[ind]     = phiij[4];
      mcs->virij[ind]    = phiij[1];
      mcs->vir12ij[ind]  = phiij[3];
      mcs->vir6ij[ind]   = phiij[5];

      mcsETrial            += mcs->eij[ind];
      mcsE12Trial            += mcs->e12ij[ind];
      mcsE6Trial            += mcs->e6ij[ind];
      mcsVirTrial          += mcs->virij[ind];
      mcsVir12Trial          += mcs->vir12ij[ind];
      mcsVir6Trial          += mcs->vir6ij[ind];
    }
    
    if ( fabs( mcs->eij[ind] - eijTest[ind] ) > 0.00001 ) {
      printf("r[%lu]: %.5G  r[%lu]: %.5G  rij[%lu][%lu]: %.5G  rijTest[%lu][%lu]: %.5G  eijTest[%lu][%lu] %.5G   eij[%lu][%lu] = %.5G\n",ii,mcs->r[ii],jj,mcs->r[jj],ii,jj,mcs->rij[ind],ii,jj,rijTest[ind],ii,jj,eijTest[ind],ii,jj,mcs->eij[ind]);
    }
  }
  #pragma omp single
  {
  mcs->E             = mcsETrial;
  mcs->E12           = mcsE12Trial;
  mcs->E6            = mcsE6Trial;
  mcs->Vir           = mcsVirTrial;
  mcs->Vir12         = mcsVir12Trial;
  mcs->Vir6          = mcsVir6Trial;
  }
  }
  
  else {
    #pragma omp single
    {
    resultFlag = 0;
    fprintf(stdout,"####### Energy was just verified ########  E = %.8G ######## ETest = %.8G #######\n",mcs->E,mcsETest);
    }
  }
  
  return resultFlag;
}


// Adjust maximum trial displacement
//   The algorithm is from Swendsen, Physics Procedia 15.
int maxDisAdjust(struct MCState *mcs) {
  static double idealRatio = 0.5;
  double actualRatio;
  
  actualRatio = (double) mcs->dAcc[0] / (mcs->dAcc[0] + mcs->dAcc[1]);
  mcs->maxStep = mcs->maxStep*log(0.672924*idealRatio + 0.0644284)/log(0.672924*(actualRatio + 0.0644284));
  if ( mcs->maxStep < 0.002 ) {
    mcs->maxStep = 0.002;
  }
  else if ( mcs->maxStep > 0.5) {
    mcs->maxStep = 0.5;
  }
  printf("Step: %lu  Updating max Step...dAcc: %lu,%lu   new maxStep: %.5G\n",mcs->sn,mcs->dAcc[0],mcs->dAcc[1],mcs->maxStep);
  
  return 0;
}


// Adjust maximum trial volume change
//   The algorithm is from Swendsen, Physics Procedia 15.
int maxDVAdjust(struct MCState *mcs) {
  static unsigned long int vAErrNtot = 0;
  static double idealRatio = 0.5;
  double actualRatio;
  
  if ( (mcs->vAcc[0] + mcs->vAcc[1] - vAErrNtot) > 0 ) {
    vAErrNtot = mcs->vAcc[0] + mcs->vAcc[1];
    actualRatio = (double) mcs->vAcc[0] / (mcs->vAcc[0] + mcs->vAcc[1]);
    mcs->maxdl = mcs->maxdl*log(0.672924*idealRatio + 0.0644284)/log(0.672924*(actualRatio + 0.0644284));
    if ( mcs->maxdl < 0.002*mcs->N) {
      mcs->maxdl = 0.002*mcs->N;
    }
    else if ( mcs->maxdl > 0.10*mcs->N) {
      mcs->maxdl = 0.50*mcs->N;
    }
    printf("Step: %lu  Updating max vol.ch....vAcc: %lu,%lu   new maxStep: %.5G\n",mcs->sn,mcs->vAcc[0],mcs->vAcc[1],mcs->maxdl);
  }
  
  return 0;
}


// Print the system's current energy to stdout
int printE(struct MCState *mcs) {
  printf("\nE = %.8G\n",mcs->E);
  return 0;
}


// Print the cumulative acceptances and rejections for each trial type
//  to stdout
int printAcc(struct MCState *mcs) {
  printf("Accepted/Rejected: Displacements VolumeChanges\n              \
          %lu/%lu          %lu/%lu\n",mcs->dAcc[0],mcs->dAcc[1],mcs->vAcc[0], \
          mcs->vAcc[1]);
  return 0;
}


// Full attempt at volume change trial (not using knowledge of previous state
//  to short-cut the calculation). Update system state according to outcome.
int fav(struct MCState *mcs) {
	unsigned long int ii,jj,ind;
	double rij1a,rij3a,rij6a,rij7a,rij12a,rij13a;
	double phiij[6];
	
	#pragma omp single
	{
	mcs->bndchk = 0;
	mcsETrial = 0;
	mcsVirTrial = 0;
	mcs->dl = (mcs->rn-0.5)*2*mcs->maxdl;
	
	params[0] = 1;
	params[1] = mcs->l + mcs->dl;
	
	lRat1 = params[1] / mcs->l;
	for (ind = 0; ind < mcs->N; ind++) {
		mcs->rTrial[ind] = mcs->r[ind]*lRat1;
	}

	mcsETrial     = 0;
	mcsE12Trial   = 0;
	mcsE6Trial    = 0;
	mcsVirTrial   = 0;
	mcsVir12Trial = 0;
	mcsVir6Trial  = 0;
	}
	
	#pragma omp for private(ii,jj,ind,rij1a,rij3a,rij6a,rij7a,rij12a,rij13a) \
             reduction (+:mcsETrial,mcsE12Trial,mcsE6Trial,mcsVirTrial,mcsVir12Trial,mcsVir6Trial)
	for (ind = 0; ind < mcs->numPairs; ind++) {
		ii            = mcs->iii[ind];
		jj            = mcs->jjj[ind];

		if ( mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
			mcs->eijTrial[ind]        = 0;
			mcs->e12ijTrial[ind]      = 0;
			mcs->e6ijTrial[ind]       = 0;
			mcs->virijTrial[ind]      = 0;
			mcs->vir12ijTrial[ind]    = 0;
			mcs->vir6ijTrial[ind]     = 0;
		}
		else {
			rij1a          = mcs->rTrial[jj] - mcs->rTrial[ii];
			
			mcs->phi(&rij1a,mcs->potCutOff,params,phiij);
			mcs->eijTrial[ind]      = phiij[0];
			mcs->e12ijTrial[ind]    = phiij[2];
			mcs->e6ijTrial[ind]     = phiij[4];
			mcs->virijTrial[ind]    = phiij[1];
			mcs->vir12ijTrial[ind]  = phiij[3];
			mcs->vir6ijTrial[ind]   = phiij[5];
	
			mcsETrial            += mcs->eijTrial[ind];
			mcsE12Trial            += mcs->e12ijTrial[ind];
			mcsE6Trial            += mcs->e6ijTrial[ind];
			mcsVirTrial          += mcs->virijTrial[ind];
			mcsVir12Trial          += mcs->vir12ijTrial[ind];
			mcsVir6Trial          += mcs->vir6ijTrial[ind];
		}
	}
	mcs->ETrial           = mcsETrial;
	mcs->E12Trial           = mcsE12Trial;
	mcs->E6Trial           = mcsE6Trial;
	mcs->VirTrial         = mcsVirTrial;
	mcs->Vir12Trial         = mcsVir12Trial;
	mcs->Vir6Trial         = mcsVir6Trial;

	#pragma omp single
	{

	mcs->bf = exp(-(mcs->ETrial - mcs->E + mcs->P*mcs->dl) / mcs->T + mcs->N*log(lRat1));

	if (mcs->bf < 1.0) {
		ran = gsl_rng_uniform(mcs->rangen);
	}

	if (mcs->bf >= 1.0 || mcs->bf > ran) {
		mcs->vAcc[0]++;
		mcs->l = mcs->l + mcs->dl;
		mcs->E = mcs->ETrial;
		mcs->E12 = mcs->E12Trial;
		mcs->E6 = mcs->E6Trial;
		mcs->Vir = mcs->VirTrial;
		mcs->Vir12 = mcs->Vir12Trial;
		mcs->Vir6 = mcs->Vir6Trial;
		for (ind = 0; ind < mcs->N; ind++) {
			mcs->r[ind] = mcs->rTrial[ind];
		}
		for (ind = 0; ind < mcs->numPairs; ind++) {
			ii = mcs->iii[ind];
			jj = mcs->jjj[ind];
			mcs->rij[ind] = mcs->r[jj] - mcs->r[ii];
			mcs->eij[ind] = mcs->eijTrial[ind];
			mcs->e12ij[ind] = mcs->e12ijTrial[ind];
			mcs->e6ij[ind] = mcs->e6ijTrial[ind];
			mcs->virij[ind] = mcs->virijTrial[ind];
			mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
			mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
		}
	}
	else {
		mcs->vAcc[1]++;
	}
	fflush(stdout);
	}
	
	return 0;
	#pragma omp barrier
}


// Quickly calculate current density and g(x), using knowledge of previous state.
int qagrho(struct MCState *mcs) {
	
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

int getRelaxFlag(struct MCState *mcs) {
    return mcs->relaxFlag;
}

bool getRestartFlag(struct MCState *mcs) {
    return mcs->isRestart;
}

// Relax volume to the minimum energy point, when considering the external
//   pressure.
int relaxVolume(struct MCState *mcs) {
//int fav(struct MCState *mcs) {
    unsigned long int ii,jj,ind;
    double relaxCrit,rij1a,rij3a,rij6a,rij7a,rij12a,rij13a;
    double phiij[6];
    int count,maxTries;
    
    maxTries = 20;
    lTryMin = 0;
    lTryMax = 1E10;
    b = 1000;

    
    #pragma omp single
    {
    printf("Relaxing the Volume...Starting Length,Energy: %.5G,%.5G\n",mcs->l,mcs->E);
    }

for ( count = 0; count < 20; count++ ) {
    // 1. Increase volume 2% (move particles to match)
    #pragma omp single
    {
    printf("Relaxation step %d...",count);
    //h = 0.02*mcs->l;
    h = 0.1;
    mcs->dl = h;
    }
    EUp = calculateEnergyOfTrialVolumeChange(mcs,mcs->dl);
   
    /*
    #pragma omp single
    {
    printf("Relaxation step %d...\n",count);
    //h = 0.02*mcs->l;
    h = 5;
    mcs->dl = h;
    
    params[0] = 0;
    
    lRat1 = (mcs->l + mcs->dl) / mcs->l;
    for (ind = 0; ind < mcs->N; ind++) {
        mcs->rTrial[ind] = mcs->r[ind]*lRat1;
    }
    
    mcsETrial     = 0;
    mcsE12Trial   = 0;
    mcsE6Trial    = 0;
    }


    // 2. Calculate "trial" energy
    #pragma omp for private(ii,jj,ind,rij1a,rij3a,rij6a,rij7a,rij12a,rij13a) \
         reduction (+:mcsETrial,mcsE12Trial,mcsE6Trial,mcsVirTrial,mcsVir12Trial,mcsVir6Trial)
    for (ind = 0; ind < mcs->numPairs; ind++) {
        ii            = mcs->iii[ind];
        jj            = mcs->jjj[ind];

        if ( mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
            mcs->eijTrial[ind]        = 0;
            mcs->e12ijTrial[ind]      = 0;
            mcs->e6ijTrial[ind]       = 0;
        }
        else {
            rij1a          = mcs->rTrial[jj] - mcs->rTrial[ii];
            
            phiij                   = mcs->phi(&rij1a,params);
            mcs->eijTrial[ind]      = phiij[0];
            mcs->e12ijTrial[ind]    = phiij[2];
            mcs->e6ijTrial[ind]     = phiij[4];
            
            mcsETrial               += mcs->eijTrial[ind];
            mcsE12Trial             += mcs->e12ijTrial[ind];
            mcsE6Trial              += mcs->e6ijTrial[ind];
        }
    }

    #pragma omp single
    {
    mcs->ETrial             = mcsETrial;
    mcs->E12Trial           = mcsE12Trial;
    mcs->E6Trial            = mcsE6Trial;
    EUp                     = mcs->ETrial;
    }
    */


// 3. Decrease volume 2% from original (move particles to match)
    #pragma omp single
    {
        //mcs->bndchk = 0;
        mcs->dl = -h;
    }
    EDown = calculateEnergyOfTrialVolumeChange(mcs,mcs->dl);

    /*    
    #pragma omp single
    {
        //mcs->bndchk = 0;
        mcs->dl = -h;
        
        params[0] = 0;
        
        lRat1 = (mcs->l + mcs->dl) / mcs->l;
        for (ind = 0; ind < mcs->N; ind++) {
            mcs->rTrial[ind] = mcs->r[ind]*lRat1;
        }

        mcsETrial     = 0;
        mcsE12Trial   = 0;
        mcsE6Trial    = 0;
        
    }

// 4. Calculate "trial" energy
	#pragma omp for private(ii,jj,ind,rij1a,rij3a,rij6a,rij7a,rij12a,rij13a) \
             reduction (+:mcsETrial,mcsE12Trial,mcsE6Trial,mcsVirTrial,mcsVir12Trial,mcsVir6Trial)
	for (ind = 0; ind < mcs->numPairs; ind++) {
		ii            = mcs->iii[ind];
		jj            = mcs->jjj[ind];

		if ( mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
			mcs->eijTrial[ind]        = 0;
			mcs->e12ijTrial[ind]      = 0;
			mcs->e6ijTrial[ind]       = 0;
		}
		else {
			rij1a          = mcs->rTrial[jj] - mcs->rTrial[ii];
			
			phiij                   = mcs->phi(&rij1a,params);
			mcs->eijTrial[ind]      = phiij[0];
			mcs->e12ijTrial[ind]    = phiij[2];
			mcs->e6ijTrial[ind]     = phiij[4];
	
			mcsETrial               += mcs->eijTrial[ind];
			mcsE12Trial             += mcs->e12ijTrial[ind];
			mcsE6Trial              += mcs->e6ijTrial[ind];
		}
	}

        #pragma omp single
        {
	mcs->ETrial             = mcsETrial;
	mcs->E12Trial           = mcsE12Trial;
	mcs->E6Trial            = mcsE6Trial;
	EDown                   = mcs->ETrial;
        */


// 5. From second-order finite difference, estimate the volume at which
//      -dE/dL = P
        #pragma omp single
        {
        double first,second;
        double relaxMax;

        first  = (EUp - EDown) / (2*h);
        second = (EUp - 2.0*mcs->E + EDown) / (h*h);
        dlEstimate = -( mcs->P - ( (double) mcs->N/mcs->l)*mcs->T  + first) / second;
	//printf("EUp,E,EDown,first,second,dlEstimate = %.5G,%.5G,%.5G,%.5G,%.5G,%.5G\n", \
        //          EUp,mcs->E,EDown,first,second,dlEstimate);
        relaxMax = (double) 0.10*mcs->N;
        if ( fabs(dlEstimate) > relaxMax ) {
            if ( dlEstimate < 0 ) {
                dlEstimate = -relaxMax;
            }
            else {
                dlEstimate = relaxMax;
            }
        }
        
        if ( mcs->l + dlEstimate > lTryMax ) {
            dlEstimate = 0.5*(lTryMax - mcs->l);
        }
        else if ( mcs->l + dlEstimate < lTryMin) {
            dlEstimate = 0.5*(lTryMin - mcs->l);
        }
        
        if ( dlEstimate > 0.0 ) {
            lTryMin = mcs->l;
        }
        else {
            lTryMax = mcs->l;
        }
        
        }
// 6. If the current volume is within relaxCrit of the estimated "minimum energy"
//      volume, move and exit the loop.  Otherwise, move and go back to
//      1.
        relaxCrit = 0.0025*mcs->N;
        if ( fabs(dlEstimate) < relaxCrit ) {
            moveVolume(mcs,mcs->l + dlEstimate);
            
            #pragma omp single
            {
	    printf("Relaxation converged. Length,energy: %.10G,%.10G\n",mcs->l,mcs->E);
            }
            return 0;
        }
	else {
            moveVolume(mcs,mcs->l + dlEstimate);
            
            #pragma omp single
            {
            printf("length,energy: %.5g,%.5G\n",mcs->l,mcs->E);
            }
      /*
   #pragma omp single
            {
         
            lRat1 = (mcs->l + dlEstimate) / mcs->l;
            mcs->l = mcs->l + dlEstimate;
            printf("Adjusting length to %.5G...",mcs->l);
            fflush(stdout);
            for (ind = 0; ind < mcs->N; ind++) {
                mcs->r[ind] = mcs->r[ind]*lRat1;
            }
            
            params[0] = 1;
            params[1] = mcs->l;
            
            mcsETrial     = 0;
            mcsE12Trial   = 0;
            mcsE6Trial    = 0;
            mcsVirTrial   = 0;
            mcsVir12Trial = 0;
            mcsVir6Trial  = 0;
            }

	    #pragma omp for private(ii,jj,ind,rij1a,rij3a,rij6a,rij7a,rij12a,rij13a) \
                reduction (+:mcsETrial,mcsE12Trial,mcsE6Trial,mcsVirTrial,mcsVir12Trial,mcsVir6Trial)
            for (ind = 0; ind < mcs->numPairs; ind++) {
		ii            = mcs->iii[ind];
		jj            = mcs->jjj[ind];

		if ( mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
			mcs->eij[ind]        = 0;
			mcs->e12ij[ind]      = 0;
			mcs->e6ij[ind]       = 0;
			mcs->virij[ind]      = 0;
			mcs->vir12ij[ind]    = 0;
			mcs->vir6ij[ind]     = 0;
		}
		else {
			mcs->rij[ind]      = mcs->r[jj] - mcs->r[ii];
			
			phiij              = mcs->phi(&(mcs->rij[ind]),params);
			mcs->eij[ind]      = phiij[0];
			mcs->e12ij[ind]    = phiij[2];
			mcs->e6ij[ind]     = phiij[4];
			mcs->virij[ind]    = phiij[1];
			mcs->vir12ij[ind]  = phiij[3];
			mcs->vir6ij[ind]   = phiij[5];
	
			mcsETrial          += mcs->eij[ind];
			mcsE12Trial        += mcs->e12ij[ind];
			mcsE6Trial         += mcs->e6ij[ind];
			mcsVirTrial        += mcs->virij[ind];
			mcsVir12Trial      += mcs->vir12ij[ind];
			mcsVir6Trial       += mcs->vir6ij[ind];
		}
            }
            
            #pragma omp single
            {
            mcs->E           = mcsETrial;
            printf("Current Energy: %.3G\n",mcs->E);
            mcs->E12         = mcsE12Trial;
            mcs->E6          = mcsE6Trial;
            mcs->Vir         = mcsVirTrial;
            mcs->Vir12       = mcsVir12Trial;
            mcs->Vir6        = mcsVir6Trial;
	    }
        */    
        }
}

    #pragma omp single
    {
    printf("WARNING: Relaxation failed to converge after 20 steps. Length,energy: "
              "%.5G,%.5G\n",mcs->l,mcs->E);
    }
    
    return 1;
    #pragma omp barrier


    
/*
// BISECTION METHOD
    

// Try increasing L
#pragma omp single
    {
    l1 = mcs->l;
    E11 = mcs->E;
    E12 = calculateEnergyOfTrialVolumeChange(mcs,h);
    dE1 = (E22-E21)/h;
    
    l3 = l1 + b;
    E31 = calculateEnergyOfTrialVolumeChange(mcs,b);
    E32 = calculateEnergyOfTrialVolumeChange(mcs,b+h);
    dE3 = (E32-E31)/h;
    
    for ( count = 0; count < maxTries; count++) {
        if ( dE2 < -mcs->P ) {
            l3 = l2 + b;
            E31 = calculateEnergyOfTrialVolumeChange(mcs,b);
            E32 = calculateEnergyOfTrialVolumeChange(mcs,b+h);
            dE3 = (E32-E31)/h;
            if ( dE3 > -mcs->P ) {
                l1 = l2;
                dE1 = dE2;

                b = 0.5*b;

                l2 = l1+b;
                moveVolume(mcs,b);
                E21 = mcs->E;
                E22 = calculateEnergyOfTrialVolumeChange(mcs,h);
                dE2 = (E22-E21)/h;
            }
            else {
                l2 = l2+b;
                moveVolume(mcs,b);
                E21 = mcs->E;
                E22 = calculateEnergyOfTrialVolumeChange(mcs,h);
                dE2 = (E22-E21)/h;
            }
        }


        d2E = (dEAbove - dEHere)/b;
        else if ( d2E > 0 )



        EDown1 = calculateEnergyofTrialVolumeChange(mcs,-b);
        EDown2 = calculateEnergyofTrialVolumeChange(mcs,-b+h);
        dEDown = (EDown2-EDown1)/h;

// 4. Calculate "trial" energy
	#pragma omp for private(ii,jj,ind,rij1a,rij3a,rij6a,rij7a,rij12a,rij13a) \
             reduction (+:mcsETrial,mcsE12Trial,mcsE6Trial,mcsVirTrial,mcsVir12Trial,mcsVir6Trial)
	for (ind = 0; ind < mcs->numPairs; ind++) {
		ii            = mcs->iii[ind];
		jj            = mcs->jjj[ind];

		if ( mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
			mcs->eijTrial[ind]        = 0;
			mcs->e12ijTrial[ind]      = 0;
			mcs->e6ijTrial[ind]       = 0;
		}
		else {
			rij1a          = mcs->rTrial[jj] - mcs->rTrial[ii];
			
			phiij                   = mcs->phi(&rij1a,params);
			mcs->eijTrial[ind]      = phiij[0];
			mcs->e12ijTrial[ind]    = phiij[2];
			mcs->e6ijTrial[ind]     = phiij[4];
	
			mcsETrial               += mcs->eijTrial[ind];
			mcsE12Trial             += mcs->e12ijTrial[ind];
			mcsE6Trial              += mcs->e6ijTrial[ind];
		}
	}

        #pragma omp single
        {
	mcs->ETrial             = mcsETrial;
	mcs->E12Trial           = mcsE12Trial;
	mcs->E6Trial            = mcsE6Trial;
	EUp                     = mcs->ETrial;
        }

        if ( mcs->E 
*/

}
    

unsigned long int getStepNum(struct MCState *mcs) {
    return mcs->sn;
}



double calculateEnergyOfTrialVolumeChange(struct MCState *mcs,double dl) {
        unsigned long int ind,ii,jj;
        double rij1a;
        double phiij[6];
        #pragma omp single
        {
        params[0] = 0;
        
        lRat1 = (mcs->l + dl) / mcs->l;
        for (ind = 0; ind < mcs->N; ind++) {
            mcs->rTrial[ind] = mcs->r[ind]*lRat1;
        }

        mcsETrial     = 0;
        mcsE12Trial   = 0;
        mcsE6Trial    = 0;
        }

	#pragma omp for private(ii,jj,ind,rij1a) \
             reduction (+:mcsETrial,mcsE12Trial,mcsE6Trial,mcsVirTrial,mcsVir12Trial,mcsVir6Trial)
	for (ind = 0; ind < mcs->numPairs; ind++) {
		ii            = mcs->iii[ind];
		jj            = mcs->jjj[ind];

		if ( mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
			mcs->eijTrial[ind]        = 0;
			mcs->e12ijTrial[ind]      = 0;
			mcs->e6ijTrial[ind]       = 0;
		}
		else {
			rij1a          = mcs->rTrial[jj] - mcs->rTrial[ii];
			
			mcs->phi(&rij1a,mcs->potCutOff,params,phiij);
			mcs->eijTrial[ind]      = phiij[0];
			mcs->e12ijTrial[ind]    = phiij[2];
			mcs->e6ijTrial[ind]     = phiij[4];
	
			mcsETrial               += mcs->eijTrial[ind];
			mcsE12Trial             += mcs->e12ijTrial[ind];
			mcsE6Trial              += mcs->e6ijTrial[ind];
		}
	}

        return mcsETrial;
}
    


int moveVolume(struct MCState *mcs,double l) {
        unsigned long int ind,ii,jj;
        double rij1a;
        double phiij[6];

	#pragma omp single
	{
	mcsETrial = 0;
	mcsVirTrial = 0;
	
	params[0] = 1;
	params[1] = l;
	
	lRat1 = params[1] / mcs->l;
	for (ind = 0; ind < mcs->N; ind++) {
		mcs->r[ind] = mcs->r[ind]*lRat1;
	}

        mcs->l        = l;
	mcsETrial     = 0;
	mcsE12Trial   = 0;
	mcsE6Trial    = 0;
	mcsVirTrial   = 0;
	mcsVir12Trial = 0;
	mcsVir6Trial  = 0;
	}
	
	#pragma omp for private(ii,jj,ind,rij1a) \
             reduction (+:mcsETrial,mcsE12Trial,mcsE6Trial,mcsVirTrial,mcsVir12Trial,mcsVir6Trial)
	for (ind = 0; ind < mcs->numPairs; ind++) {
		ii            = mcs->iii[ind];
		jj            = mcs->jjj[ind];

		if ( mcs->nbn >= 0 && (unsigned int) jj-ii > mcs->nbn) {
			mcs->eij[ind]        = 0;
			mcs->e12ij[ind]      = 0;
			mcs->e6ij[ind]       = 0;
			mcs->virij[ind]      = 0;
			mcs->vir12ij[ind]    = 0;
			mcs->vir6ij[ind]     = 0;
		}
		else {
			mcs->rij[ind]        = mcs->r[jj] - mcs->r[ii];
			
			mcs->phi(&(mcs->rij[ind]),mcs->potCutOff,params,phiij);
			mcs->eij[ind]      = phiij[0];
			mcs->e12ij[ind]    = phiij[2];
			mcs->e6ij[ind]     = phiij[4];
			mcs->virij[ind]    = phiij[1];
			mcs->vir12ij[ind]  = phiij[3];
			mcs->vir6ij[ind]   = phiij[5];
	
			mcsETrial            += mcs->eij[ind];
			mcsE12Trial            += mcs->e12ij[ind];
			mcsE6Trial            += mcs->e6ij[ind];
			mcsVirTrial          += mcs->virij[ind];
			mcsVir12Trial          += mcs->vir12ij[ind];
			mcsVir6Trial          += mcs->vir6ij[ind];
		}
	}
	mcs->E             = mcsETrial;
	mcs->E12           = mcsE12Trial;
	mcs->E6            = mcsE6Trial;
	mcs->Vir           = mcsVirTrial;
	mcs->Vir12         = mcsVir12Trial;
	mcs->Vir6          = mcsVir6Trial;

        return mcs->l;
}




