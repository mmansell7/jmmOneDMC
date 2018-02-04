
struct MCState;
struct MCState * setupMCS();
void mcsPrintStep(struct MCState *);
void mcs_fad(struct MCState *,unsigned long int *, double *);
int printMCP(struct MCState *mcs1);
int mcs_printCoords(struct MCState *mcs);
int mcs_printRho(struct MCState *mcs);
int mcs_printG(struct MCState *mcs);
int mcs_fgrho(struct MCState *mcs);
int mcs_ugrho(struct MCState *mcs);
int mcs_qad(struct MCState *mcs,unsigned long int nm,double rndm);
int mcs_qav(struct MCState *mcs,unsigned long int nm,double rndm);

struct MCInput {
  
  unsigned long int N;
  double P;
  double T;
  unsigned long int ns;
  unsigned long int cpi;
  unsigned long int tpi;
  unsigned long int gpi;
  unsigned long int rhopi;
  unsigned long int gnb;
  unsigned long int rhonb;
  unsigned long int seed;
  double maxStep;
  double maxdl;
  double * (*phi)(double *rij, void *params);
  char potStr[80];
  double rbw;
  double gsw;
  double gbw;
  int gns;
  int nbn;
};


