#include <stdio.h>
#include "readInput.h"


//struct MCInput readInput2(int argc, char *argv[],unsigned long int *N,double *P,double *T,unsigned long int *ns,double * (*phi)(double*,void*),double *ms,
//		double *mdl,unsigned long int *cpi,unsigned long int *tpi) {
struct MCInput readInput2() {

	/*	Read input arguments:
			1.  Number of Particles
			2.  Pressure
			3.  Temperature
			4.  Neighbor number limit
			5.  Number of trial steps
			6.  Potential type
			7.  Maximum particle displacement
			8.  Maximum box length change
			9.  Configuration print interval
			10. Thermo print interval
			11. Rho bin width
			12. Number of bins for rho
			13. Rho print interval
			14. g(x) segment width
			15. Number of segments for g(x)
			16. g(x) bin width
			17. Number of bins for g(x)
			18. g(x) print interval
			19. Seed (suggested use: date +'%s')

		Identified future extensions:

	*/
	
    struct MCInput inp;
    FILE *infile = fopen("INPUT","r");
    char line[100];
    const char s[8] = "[ \t\n]";
    char *token;
    unsigned long int N,nbn,numSteps,cpi,tpi,rnb,rpi,gns,gnb,gpi,seed;
    double P,T,maxStep,maxdl,rbw,gsw,gbw;
    void * pot;
    
    
    if (infile==NULL)
    {
        printf("FATAL ERROR: INPUT not found.");
    }
   
    while ( fgets(line,100,infile) != NULL) {
        // get the property identifier token
        token = strtok(line, s);
        
        if ( strncmp(token,"N",10) == 0 ) {
            // Get the first value token
            token = strtok(NULL, s);
            /* walk through remaining value tokens */
            while( token != NULL ) {
                // strtok(NULL,delim) continues tokenizing where it last left off
	        inp.N = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("N = %lu ",inp.N);
            printf("\n");
        }
        
        else if ( strncmp(token,"P",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.P = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
            printf("P = %.5G ",inp.P);
            printf("\n");
        }

        else if ( strncmp(token,"T",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.T = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
            printf("T = %.5G ",inp.T);
            printf("\n");
        }

        else if ( strncmp(token,"NBN",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.nbn = (int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("NBN = %u ",inp.nbn);
            printf("\n");
        }
        
        else if ( strncmp(token,"NUMSTEPS",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.ns = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("NUMSTEPS = %lu ",inp.ns);
            printf("\n");
        }
        
        else if ( strncmp(token,"POT",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                strncpy(inp.potStr,token,80);
                if ( strncmp(inp.potStr,"LJ",10) != 0 ) {
                    printf("Unknown potential type requested (%s). Using Lennard-Jones potential.\n",token);
                    strncpy(inp.potStr,"LJ",80);
                }
                token = strtok(NULL, s);
            }
            printf("POT = %s ",inp.potStr);
            printf("\n");
        }
        
        else if ( strncmp(token,"MAXSTEP",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.maxStep = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
            printf("MAXSTEP = %.5G ",inp.maxStep);
            printf("\n");
        }
        
        else if ( strncmp(token,"MAXDV",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.maxdl = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
            printf("MAXDV = %.5G ",inp.maxdl);
            printf("\n");
        }
       
        else if ( strncmp(token,"CPI",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.cpi = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("CPI = %lu ",inp.cpi);
            printf("\n");
        }
       
        else if ( strncmp(token,"TPI",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.tpi = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("TPI = %lu ",inp.tpi);
            printf("\n");
        }
       
        else if ( strncmp(token,"GPI",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gpi = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("GPI = %lu ",inp.gpi);
            printf("\n");
        }
       
        else if ( strncmp(token,"RHOPI",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.rhopi = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("RHOPI = %lu ",inp.rhopi);
            printf("\n");
        }
       
        else if ( strncmp(token,"RBW",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.rbw = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
            printf("RBW = %.5G ",inp.rbw);
            printf("\n");
        }
       
        else if ( strncmp(token,"RHONB",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.rhonb = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("RHONB = %lu ",inp.rhonb);
            printf("\n");
        }
       
        else if ( strncmp(token,"GBW",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gbw = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
            printf("gBW = %.5G ",inp.gbw);
            printf("\n");
        }
       
        else if ( strncmp(token,"GNB",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gnb = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("GNB = %lu ",inp.gnb);
            printf("\n");
        }
       
        else if ( strncmp(token,"GSW",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gsw = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
            printf("GSW = %.5G ",inp.gsw);
            printf("\n");
        }
       
        else if ( strncmp(token,"GNS",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gns = (int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("GNS = %d ",inp.gns);
            printf("\n");
        }
        
        else if ( strncmp(token,"SEED",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.seed = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
            printf("SEED = %lu ",inp.seed);
            printf("\n");
        }
        
        else {
            token = strtok(NULL,s);
            printf("Property command %s not understood.\n",token);
        }

    }

	/* 
	*N       = (unsigned long int) strtol(argv[1],NULL,10);
	*P       = strtod(argv[2],NULL);
	*T       = strtod(argv[3],NULL);
	nbn      = (int) strtol(argv[4],NULL,10); 
	*ns      = (unsigned long int) strtol(argv[5],NULL,10);
	if (strcmp(argv[6],"lj") == 0) {
		printf("Using Lennard-Jones potential\n");
		phi = phiLJ;
	}
	else {
		printf("Unknown potential type requested (%s). Using Lennard-Jones potential.\n",argv[6]);
		phi = phiLJ;
	}

	*ms      = (double) strtod(argv[7],NULL);
	*mdl     = (double) strtod(argv[8],NULL);
	*cpi     = (unsigned long int) strtol(argv[9],NULL,10);
	*tpi     = (unsigned long int) strtol(argv[10],NULL,10);

	rbw      = (double) strtod(argv[11],NULL);
	rhonb    = (unsigned long int) strtol(argv[12],NULL,10);
	rhopi    = (unsigned long int) strtol(argv[13],NULL,10);

        
	gsw      = (double) strtod(argv[14],NULL);
	gns      = (int) strtol(argv[15],NULL,10);
	gbw      = (double) strtod(argv[16],NULL);
	gnb      = (unsigned long int) strtol(argv[17],NULL,10);
	gpi      = (unsigned long int) strtol(argv[18],NULL,10);
	seed     = (unsigned long int) strtol(argv[19],NULL,10);

	inp.gbw = gbw;
	inp.gnb = gnb;
	inp.seed = seed;
	*/
		
        return inp;
	


}
