#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jmmMCState.h"
#include "readInput.h"


struct MCInput readInput(char *fstr) {

	/*	Read input arguments:
                        0.  Is this a restart?
			1.  Number of Particles
			2.  Pressure
			3.  Temperature
			4.  Neighbor number limit
			5.  Number of trial steps
			6.  Potential type
			7.  Initial maximum particle displacement
			8.  Initial maximum box length change
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
	
    struct MCInput inp = {.isRestart = false};
    FILE *infile = fopen(fstr,"r");
    char line[100];
    const char s[8] = "[ \t\n]";
    const int maxNumTokens = 20;
    char *tokens[maxNumTokens];
    int ii, counter;
    
    if (infile==NULL)
    {
        printf("FATAL ERROR: INPUT not found.");
    }
   
    while ( fgets(line,200,infile) != NULL) {
        // get the first word from line into tokens[0]
        
        counter = 0;
        tokens[counter] = strtok(line,s);
        while( tokens[counter] != NULL ) {
            counter++;
            tokens[counter] = strtok(NULL,s);
        }
        
        if ( strncmp(tokens[0],"RESTART",10) == 0 ) {
            const int caseStrNum = 3;
            const char *falstr[] = {"FALSE","False","false"}; // List of possible potential keywords
            inp.isRestart = true;
            for ( ii = 0; ii < caseStrNum; ii++ ) {
                if ( strncmp(tokens[1],falstr[ii],20) == 0 ) {
                    inp.isRestart = false;
                }
            }
        }
        
        else if ( strncmp(tokens[0],"N",10) == 0 ) {
	    inp.N = (unsigned long int) strtol(tokens[1],NULL,10);
        }
        
        else if ( strncmp(tokens[0],"P",10) == 0 ) {
            inp.P = (double) strtod(tokens[1],NULL);
        }

        else if ( strncmp(tokens[0],"T",10) == 0 ) {
            inp.T = (double) strtod(tokens[1],NULL);
        }

        else if ( strncmp(tokens[0],"NBN",10) == 0 ) {
            inp.nbn = (int) strtol(tokens[1],NULL,10);
        }
        
        else if ( strncmp(tokens[0],"NUMSTEPS",10) == 0 ) {
            inp.ns = (unsigned long int) strtol(tokens[1],NULL,10);
        }
        
        else if ( strncmp(tokens[0],"POT",10) == 0 ) {
            const int potLen = 3;
            const char *pots[] = {"LJcut","LJ","HARMONIC"}; // List of possible potential keywords
            int potFlag;
            
            strncpy(inp.potStr,tokens[1],80); // Copy up to 80 characters from the tokens[1] to inp.potStr
              
            // Loop through pots, confirming that the requested potential keyword is among the allowed
            //   values, or handling the situation if the keyword is not among the accepted values.
            potFlag = -1;
            for (ii = 0; ii < potLen; ii++) {
                if ( strncmp(inp.potStr,pots[ii],20) == 0 ) {
                    potFlag = ii;
                    break;
                }
            }
            if ( potFlag < 0 ) {
                printf("Unknown potential type requested (%s). Using Lennard-Jones potential with no cut-off.\n",tokens[1]);
                strncpy(inp.potStr,"LJ",80);
            }
            
            inp.potCutOff = strtod(tokens[2],NULL);
        }
        
        else if ( strncmp(tokens[0],"MAXSTEP",10) == 0 ) {
            inp.maxStep = (double) strtod(tokens[1],NULL);
        }
        
        else if ( strncmp(tokens[0],"MAXDV",10) == 0 ) {
            inp.maxdl = (double) strtod(tokens[1],NULL);
        }
       
        else if ( strncmp(tokens[0],"CPI",10) == 0 ) {
            inp.cpi = (unsigned long int) strtol(tokens[1],NULL,10);
        }
       
        else if ( strncmp(tokens[0],"TPI",10) == 0 ) {
            inp.tpi = (unsigned long int) strtol(tokens[1],NULL,10);
        }
       
        else if ( strncmp(tokens[0],"GPI",10) == 0 ) {
            inp.gpi = (unsigned long int) strtol(tokens[1],NULL,10);
        }
       
        else if ( strncmp(tokens[0],"RHOPI",10) == 0 ) {
            inp.rhopi = (unsigned long int) strtol(tokens[1],NULL,10);
        }
       
        else if ( strncmp(tokens[0],"RBW",10) == 0 ) {
            inp.rbw = (double) strtod(tokens[1],NULL);
        }
       
        else if ( strncmp(tokens[0],"RHONB",10) == 0 ) {
            inp.rhonb = (unsigned long int) strtol(tokens[1],NULL,10);
        }
       
        else if ( strncmp(tokens[0],"GBW",10) == 0 ) {
            inp.gbw = (double) strtod(tokens[1],NULL);
        }
       
        else if ( strncmp(tokens[0],"GNB",10) == 0 ) {
            inp.gnb = (unsigned long int) strtol(tokens[1],NULL,10);
        }
       
        else if ( strncmp(tokens[0],"GSW",10) == 0 ) {
            inp.gsw = (double) strtod(tokens[1],NULL);
        }
       
        else if ( strncmp(tokens[0],"GNS",10) == 0 ) {
            inp.gns = (int) strtol(tokens[1],NULL,10);
        }
        
        else if ( strncmp(tokens[0],"SEED",10) == 0 ) {
            inp.seed = (unsigned long int) strtol(tokens[1],NULL,10);
        }
        
        else if ( strncmp(tokens[0],"ENGCHECK",10) == 0 ) {
            inp.eci = (unsigned long int) strtol(tokens[1],NULL,10);
        }
        
        else if ( strncmp(tokens[0],"DADJ",10) == 0 ) {
            inp.mdai = (unsigned long int) strtol(tokens[1],NULL,10);
        }
        
        else if ( strncmp(tokens[0],"VADJ",10) == 0 ) {
            inp.mvai = (unsigned long int) strtol(tokens[1],NULL,10);
        }
        
        else if ( strncmp(tokens[0],"RELAX",10) == 0 ) {
            inp.relaxFlag = 1;
        }
        
        else {
            printf("Property command %s not understood.\n",tokens[0]);
        }

    }
    
    return inp;

}


int printInput(struct MCInput inp) {

    printf("N = %lu\n",inp.N);
    printf("P = %.5G\n",inp.P);
    printf("T = %.5G\n",inp.T);
    if ( inp.nbn > 0 ) {
        printf("NBN = %d\n",inp.nbn);
    }
    else {
        printf("NBN = No neighbor number limit.\n");
    }
    printf("NUMSTEPS = %lu\n",inp.ns);
    printf("POT = %s, CUTOFF = %.5G\n",inp.potStr,inp.potCutOff);
    printf("MAXSTEP = %.5G\n",inp.maxStep);
    printf("MAXDV = %.5G\n",inp.maxdl);
    printf("ENGCHECK = %lu\n",inp.eci);
    printf("DADJ = %lu\n",inp.mdai);
    printf("VADJ = %lu\n",inp.mvai);
    printf("CPI = %lu\n",inp.cpi);
    printf("TPI = %lu\n",inp.tpi);
    printf("GPI = %lu\n",inp.gpi);
    printf("RHOPI = %lu\n",inp.rhopi);
    printf("RBW = %.5G\n",inp.rbw);
    printf("RHONB = %lu\n",inp.rhonb);
    printf("GBW = %.5G\n",inp.gbw);
    printf("GNB = %lu\n",inp.gnb);
    printf("GSW = %.5G\n",inp.gsw);
    printf("GNS = %d\n",inp.gns);
    printf("SEED = %lu\n",inp.seed);
    printf("RELAXATION FLAG = %d\n",inp.relaxFlag);

    return 0;
}
