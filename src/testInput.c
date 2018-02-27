#include <stdio.h>
#include "jmmMCState.h"
#include "readInput.h"

int main(int argc, char **argv[]) {
  printf("test 1: %d\n",test1());
  return 0;
}

int test1() {

  char *fstr1 = "TEST1INPUT";
  FILE *f = fopen(fstr1,"w");

  unsigned long int N        = 5;
  double P                   = 0.5;
  double T                   = 0.5;
  unsigned long int NUMSTEPS = 2;
  char *POT                  = "YT";
  int NBN                    = -1;
  double MAXSTEP             = 0.1;
  double MAXDV               = 1.0;
  unsigned long int CPI      = 1;
  unsigned long int TPI      = 1;
  double RBW                 = 0.01;
  unsigned long int RHONB    = 1000;
  unsigned long int RHOPI    = 1;
  double GSW                 = 1;
  unsigned int GNS           = 5;
  double GBW                 = 0.01;
  unsigned long int GNB      = 1000;
  unsigned long int GPI      = 1;
  unsigned long int SEED     = 125;
  double YADA                = 16.543;
  unsigned long int ENGCHECK = 1;
  unsigned long int DADJ     = 1;
  unsigned long int VADJ     = 1;

  fprintf(f,"N      %lu\n"
          "P      %.5G\n"
          "T   %.5G\n"
          "NUMSTEPS\t\t%lu\n"
          "POT\t    %s\n"
          "NBN\t\t%d\n"
          "MAXSTEP    %.5G\n"
          "MAXDV     %.5G\n"
          "CPI\t\t   %lu\n"
          "TPI\t\t%lu\n"
          "RBW\t\t%.5G\n"
          "RHONB\t\t%lu\n"
          "RHOPI\t\t%lu\n"
          "GSW\t\t%.5G\n"
          "GNS\t\t%u\n"
          "GBW\t\t%.5G\n"
          "GNB\t\t%lu\n"
          "GPI\t\t%lu\n"
          "SEED\t\t%lu\n"
          "YADA\t\t%.5G\n"
          "ENGCHECK\t\t%lu\n"
          "DADJ\t\t%lu\n"
          "VADJ\t\t%lu",
          N,P,T,NUMSTEPS,POT,NBN,MAXSTEP,MAXDV,CPI,TPI,RBW,RHONB,RHOPI,
          GSW,GNS,GBW,GNB,GPI,SEED,YADA,ENGCHECK,DADJ,VADJ);

  fclose(f);
  
  struct MCInput mci1 = readInput(fstr1);
  
  int Ncheck = !(mci1.N == N);
  int Pcheck = !(mci1.P == P);
  int Tcheck = !(mci1.T == T);
  int NUMSTEPScheck = !(mci1.ns == NUMSTEPS);
  int POTcheck = !(strncmp(mci1.potStr,POT,10));
  int NBNcheck = !(mci1.nbn == NBN);
  int MAXSTEPcheck = !(mci1.maxStep == MAXSTEP);
  int MAXDVcheck = !(mci1.maxdl == MAXDV);
  int CPIcheck = !(mci1.cpi == CPI);
  int TPIcheck = !(mci1.tpi == TPI);
  int RBWcheck = !(mci1.rbw == RBW);
  int RHONBcheck = !(mci1.rhonb == RHONB);
  int RHOPIcheck = !(mci1.rhopi == RHOPI);
  int GSWcheck = !(mci1.gsw == GSW);
  int GNScheck = !(mci1.gns == GNS);
  int GBWcheck = !(mci1.gbw == GBW);
  int GNBcheck = !(mci1.gnb == GNB);
  int GPIcheck = !(mci1.gpi == GPI);
  int SEEDcheck = !(mci1.seed == SEED);
  int ENGCHECKcheck = !(mci1.eci == ENGCHECK);
  int DADJcheck = !(mci1.mdai == DADJ);
  int VADJcheck = !(mci1.mvai == VADJ);

  int check = Ncheck + 2*Pcheck + 4*Tcheck + 8*NUMSTEPScheck + 16*POTcheck + 32*NBNcheck + 64*MAXSTEPcheck + 
          128*MAXDVcheck + 256*CPIcheck + 512*TPIcheck + 1024*RBWcheck + 2048*RHONBcheck + 4096*RHOPIcheck + 
          8192*GSWcheck + 16384*GNScheck + 2*16384*GBWcheck + 4*16384*GNBcheck + 8*16384*GPIcheck + 
          16*16384*SEEDcheck + 32*16384*ENGCHECKcheck + 64*16384*DADJcheck + 128*16384*VADJcheck;

  return check;
}
