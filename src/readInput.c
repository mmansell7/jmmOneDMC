#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jmmMCState.h"
#include "readInput.h"


struct MCInput readInput() {

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
        }
        
        else if ( strncmp(token,"P",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.P = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
        }

        else if ( strncmp(token,"T",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.T = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
        }

        else if ( strncmp(token,"NBN",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.nbn = (int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
        
        else if ( strncmp(token,"NUMSTEPS",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.ns = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
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
        }
        
        else if ( strncmp(token,"MAXSTEP",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.maxStep = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
        }
        
        else if ( strncmp(token,"MAXDV",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.maxdl = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"CPI",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.cpi = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"TPI",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.tpi = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"GPI",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gpi = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"RHOPI",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.rhopi = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"RBW",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.rbw = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"RHONB",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.rhonb = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"GBW",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gbw = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"GNB",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gnb = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"GSW",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gsw = (double) strtod(token,NULL);
                token = strtok(NULL, s);
            }
        }
       
        else if ( strncmp(token,"GNS",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.gns = (int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
        
        else if ( strncmp(token,"SEED",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.seed = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
        
        else if ( strncmp(token,"ENGCHECK",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.eci = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
        
        else if ( strncmp(token,"DADJ",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.mdai = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
        
        else if ( strncmp(token,"VADJ",10) == 0 ) {
            token = strtok(NULL, s);
            while( token != NULL ) {
                inp.mvai = (unsigned long int) strtol(token,NULL,10);
                token = strtok(NULL, s);
            }
        }
        
        else {
            printf("Property command %s not understood.\n",token);
        }

    }
    
    /*
    printf("N = %lu\n",inp.N);
    printf("P = %.5G\n",inp.P);
    printf("T = %.5G\n",inp.T);
    printf("NBN = %u\n",inp.nbn);
    printf("NUMSTEPS = %lu\n",inp.ns);
    printf("POT = %s\n",inp.potStr);
    printf("MAXSTEP = %.5G\n",inp.maxStep);
    printf("MAXDV = %.5G\n",inp.maxdl);
    printf("ENGCHECK = %.5G\n",inp.eci);
    printf("DADJ = %.5G\n",inp.mdai);
    printf("VADJ = %.5G\n",inp.mvai);
    printf("CPI = %lu\n",inp.cpi);
    printf("TPI = %lu\n",inp.tpi);
    printf("GPI = %lu\n",inp.gpi);
    printf("RHOPI = %lu\n",inp.rhopi);
    printf("RBW = %.5G\n",inp.rbw);
    printf("RHONB = %lu\n",inp.rhonb);
    printf("gBW = %.5G\n",inp.gbw);
    printf("GNB = %lu\n",inp.gnb);
    printf("GSW = %.5G\n",inp.gsw);
    printf("GNS = %d\n",inp.gns);
    printf("SEED = %lu\n",inp.seed);
    */
    
    return inp;
	


}
