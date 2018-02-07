
struct MCState;
struct MCState * setupMCS();
void printStep(struct MCState *);
int printMCP(struct MCState *mcs1);
int printThermo(struct MCState *mcs);
int printCoords(struct MCState *mcs);
int printRho(struct MCState *mcs);
int printG(struct MCState *mcs);
int fgrho(struct MCState *mcs);
int ugrho(struct MCState *mcs);
int qagrho(struct MCState *mcs);
int NPTtrial(struct MCState *mcs);
int fad(struct MCState *,unsigned long int *, double *);
int qad(struct MCState *mcs,unsigned long int *nm, double *d);
int qad2(struct MCState *mcs,unsigned long int *nm, double *d);
int fav(struct MCState *mcs);
int qav(struct MCState *mcs);
int qavLJ(struct MCState *mcs);
unsigned long int incrementStep(struct MCState *mcs);
int Step(struct MCState *mcs);
int isCoordPrint(struct MCState *mcs);
int isThermoPrint(struct MCState *mcs);
int isRhoPrint(struct MCState *mcs);
int isGPrint(struct MCState *mcs);
int isECheck(struct MCState *mcs);
int isMaxDisAdjust(struct MCState *mcs);
int isMaxDVAdjust(struct MCState *mcs);
int updateThermo(struct MCState *mcs);
int ECheck(struct MCState *mcs);
int MaxDisAdjust(struct MCState *mcs);
int MaxDVAdjust(struct MCState *mcs);

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
  char potStr[80];
  double rbw;
  double gsw;
  double gbw;
  int gns;
  int nbn;
};




