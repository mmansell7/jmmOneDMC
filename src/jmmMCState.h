
struct MCState;
struct MCState * setupMCS();
void mcsPrintStep(struct MCState *);
int printMCP(struct MCState *mcs1);
int mcs_printThermo(struct MCState *mcs);
int mcs_printCoords(struct MCState *mcs);
int mcs_printRho(struct MCState *mcs);
int mcs_printG(struct MCState *mcs);
int mcs_fgrho(struct MCState *mcs);
int mcs_ugrho(struct MCState *mcs);
int mcs_qagrho(struct MCState *mcs);
int mcs_NPTtrial(struct MCState *mcs);
void mcs_fad(struct MCState *,unsigned long int *, double *);
int mcs_qad(struct MCState *mcs);
int mcs_fav(struct MCState *mcs);
int mcs_qav(struct MCState *mcs);
unsigned long int mcs_incrementStep(struct MCState *mcs);
int mcs_Step(struct MCState *mcs);
int mcs_isCoordPrint(struct MCState *mcs);
int mcs_isThermoPrint(struct MCState *mcs);
int mcs_isRhoPrint(struct MCState *mcs);
int mcs_isGPrint(struct MCState *mcs);
int mcs_isECheck(struct MCState *mcs);
int mcs_isMaxDisAdjust(struct MCState *mcs);
int mcs_isMaxDVAdjust(struct MCState *mcs);
int mcs_updateThermo(struct MCState *mcs);
int mcs_ECheck(struct MCState *mcs);
int mcs_MaxDisAdjust(struct MCState *mcs);
int mcs_MaxDVAdjust(struct MCState *mcs);

struct MCInput {
  
  unsigned long int N;
  double P;
  double T;
  unsigned long int ns;
  unsigned long int cpi;
  unsigned long int tpi;
  unsigned long int eci;
  unsigned long int mdai;
  unsigned long int mvai;
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




