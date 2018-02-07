/*
 ============================================================================
 Name        : Main.c
 Author      : Matt Mansell, NC State University, Dept. of Chemical &
                 Biomolecular Engineering
 Version     : 0.0.1
 Copyright   : Your copyright notice
 Description : Main function for One-Dimensional NPT Simulations
 ============================================================================
 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "jmmMCState.h"
#include "readInput.h"
#include "pot.h"

char *verStr = "1.1.0";
char *verDateStr = "02/05/2018";



int main (int argc, char *argv[]) {
    
    printf("########################################\n");
    printf("#       Matt Mansell's MC Code         #\n");
    printf("#      Version %s : %s      #\n",verStr,verDateStr);
    printf("########################################\n\n");
    
    // Read input from file INPUT
    struct MCInput inp = readInput("INPUT");
    
    // Using values read into inp, create an MCState struct that defines
    //  the state of the simulation throughout the simulation. But
    //  this struct is accessible through the pointer mcs, so that it
    //  can be passed by value (not by reference) to various functions
    //  that modify it - and hence, modify that state of the system
    //  under simulation.
    struct MCState *mcs = setupMCS(inp);
    
    fprintf(stdout,"Setup completed\n");
    
// Execute everything within the following block in parallel, sharing the mcs
//  struct pointer among all threads.
#pragma omp parallel default (shared) shared(mcs)
{
    // Execute the following block with only a single thread.
    #pragma omp single
    {
    printf("Running with %d threads...\n\n",omp_get_num_threads());
    
    mcsPrintStep(mcs);
    fflush(stdout);
    }
    unsigned long int ddd = 0;
    double nnn = 0.5;
    mcs_fad(mcs,&ddd,&nnn);
    
    #pragma omp sections
    {
        #pragma omp section
        {
        mcs_printCoords(mcs);
        }
        
        #pragma omp section
        {
        mcs_printRho(mcs);
        }
        
        #pragma omp section
        {
        mcs_updateThermo(mcs);
        mcs_printThermo(mcs);
        }
    }
    
    mcs_printG(mcs);
    
    // ############################### MC LOOP #######################################
    while ( mcs_incrementStep(mcs) ) {
        #pragma omp barrier
        mcs_Step(mcs);
        #pragma omp barrier
        printf("Out of mcs_Step...\n");
        #pragma omp barrier
        fflush(stdout);
        
        //#pragma omp barrier
        //#pragma omp sections private(ii,jj,dj,ind) nowait
        {
        //#pragma omp section
        #pragma omp single
        {
        if ( mcs_isCoordPrint(mcs) ) {
            mcs_printCoords(mcs);
        }
        }

        //#pragma omp section
        #pragma omp single
        {
        if ( mcs_isThermoPrint(mcs) ) {
            mcs_printThermo(mcs);
        }
        }
        
        //#pragma omp section
        #pragma omp single
        {
        if ( mcs_isRhoPrint(mcs) ) {
            mcs_printRho(mcs);
        }
        }
        
        //#pragma omp section
        #pragma omp single
        {
        if ( mcs_isMaxDisAdjust(mcs) ) {
             mcs_maxDisAdjust(mcs);
        }
        }
        
        //#pragma omp section
        #pragma omp single
        {
        if ( mcs_isMaxDVAdjust(mcs) ) {
            mcs_maxDVAdjust(mcs);
        }
        }
        }
        
        if ( mcs_isGPrint(mcs) ) {
            mcs_printG(mcs);
        }
        
        fflush(stdout);
        
    }     //  ############################### END MC LOOP #######################################

#pragma omp barrier

}

	printf("\nPROGRAM COMPLETED SUCCESSFULLY!\n");
	mcs_printE(mcs);
        mcs_printAcc(mcs);
        return 0;

}





