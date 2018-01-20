
struct MCState;
struct MCState * setupMCS();
void mcsPrintStep(struct MCState *);


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
  double rbw;
  double gsw;
  double gbw;
  int gns;
  
};


