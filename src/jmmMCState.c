#include <stdio.h>
#include "jmmMCState.h"

struct jmmMCState {
  
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

  double P;                        // Pressure (reduced units)
  double T;                        // Temperature (reduced units)
  double maxStep;                  // Maximum trial displacement length
  double maxdl;                    // Maximum volume change trial
  double * (*phi)(double);         // Pointer to the potential function
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
  
};



