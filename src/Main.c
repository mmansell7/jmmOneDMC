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
    
    // Print the step number to stdout
    printStep(mcs);
    }
    
    // In order to calculate initial energy, etc., step 0 is a "trial move"
    //  in which mcs_fad is called with arguments indicating particle 0
    //  should be moved by an amount equivalent to 50% probability within
    //  the trial displacement probability distribution (which is equal to
    //  a trial displacement of precisely 0).
    unsigned long int nnn = 0;
    double ddd = 0.5;
    fad(mcs,&nnn,&ddd);
    
    // Volume is relaxed, if this was requested in the input file.
    if ( getRelaxFlag(mcs) > 0 ) {
        relaxVolume(mcs);
    }
    
    // The following "sections" are each executed by a single thread, in 
    //  parallel with the other sections.
    #pragma omp sections
    {
        #pragma omp section
        {
        // Print coordinates at current step to the configuration file.
        printCoords(mcs);
        }
        
        #pragma omp section
        {
        // Print density bin populations at current step to the rho file.
        printRho(mcs);
        }
        
        #pragma omp section
        {
        // Accumulate thermodynamic parameters at step 0.  Accumulators are used to
        //   calculate block averages at the end of each block. 
        updateThermo(mcs);

        // Print step 0 thermodynamics to the thermo output file.
        printThermo(mcs);
        }
    }
    
    // Print g(r) or two-particle density bin populations to the g#.dat files.
    //   Note that the mcs_printG function is executed by all threads in parallel.
    printG(mcs);
    
    
    // ############################### MC LOOP #######################################
    
           // In mcs_incrementStep, a single thread attempts to increment the step
           //  number, returning 1 if the maximum step number has been
           //  reached, and 0 otherwise.  All other threads wait until the check
           //  has completed, so that all threads always return the same value.
    while ( incrementStep(mcs) ) {
        #pragma omp barrier
        // Perform a Monte Carlo step, updating the systems state appropriately.
        Step(mcs);
        #pragma omp barrier
        
        //#pragma omp barrier
        
        //In the following sections, output files are updated, and trial maximum
        //  values (displacement, volume change, etc.) are updated at the steps
        //  defined in the input file.
        //#pragma omp sections private(ii,jj,dj,ind) nowait
        {
        //#pragma omp section
        #pragma omp single
        {
        if ( isCoordPrint(mcs) ) {
            printCoords(mcs);
        }
        }

        //#pragma omp section
        #pragma omp single
        {
        if ( isThermoPrint(mcs) ) {
            printThermo(mcs);
        }
        }
        
        //#pragma omp section
        #pragma omp single
        {
        if ( isRhoPrint(mcs) ) {
            printRho(mcs);
        }
        }
        
        //#pragma omp section
        #pragma omp single
        {
        if ( isMaxDisAdjust(mcs) ) {
             maxDisAdjust(mcs);
        }
        }
        
        //#pragma omp section
        #pragma omp single
        {
        if ( isMaxDVAdjust(mcs) ) {
            maxDVAdjust(mcs);
        }
        }
        }
        
        if ( isGPrint(mcs) ) {
            printG(mcs);
        }
        
        // Volume is relaxed, if this was requested in the input file.
        if ( (getStepNum(mcs) % 10000 == 0) && (getStepNum(mcs) < 1E6) && \
               (getRelaxFlag(mcs) > 0) ) {
            relaxVolume(mcs);
        }
    
        fflush(stdout);
        
    }     //  ############################### END MC LOOP #######################################

#pragma omp barrier

}

	printf("\nPROGRAM COMPLETED SUCCESSFULLY!\n");
	printE(mcs);
        printAcc(mcs);
	freeMCS(mcs);
        return 0;

}





