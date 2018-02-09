#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "jmmMCState.h"

// Pair energy is always returned as the first element of
//  phi, and pair configurational pressure is always 
//  returned as the second.
//
// Generally, param[0] determines whether the pressure
//  contribution is calculated (0: no, 1: yes), and param[1]
//  gives the box volume, which is necessary to compute the
//  pressure contributions.


double * phiLJ(double *d, void *params) {
	double rij1,rij3,rij6,rij7,rij12,rij13;
	double phitot,phi6,phi12;
        double l;
	double virtot,vir6,vir12;
	double *phi = (double *) malloc(6*sizeof(double));
	
        double *p = (double *) params;
        
	rij1 = *d;
	rij3 = rij1*rij1*rij1;
	rij6 = 1/(rij3*rij3);
	rij12 = rij6*rij6;

	phi6 = 4*rij6;
	phi12 = 4*rij12;
        phitot = phi12 - phi6;
	phi[0] = phitot;
	phi[2] = phi12;
        phi[4] = phi6;

        if ( p[0] > 0.1 ) {
	  l      = p[1];
          rij7   = rij6/rij1;
	  rij13  = rij12/rij1;
	  
          vir6   = 24/l*rij6;
	  vir12  = 48/l*rij12;
	  virtot = vir12 - vir6;
          phi[1] = virtot;
	  phi[3] = vir12;
          phi[5] = vir6;
	}
	else {
	  phi[1] = 0;
	  phi[3] = 0;
	  phi[5] = 0;
	}

	return phi;
}



double * phiHarmonic(double *d, void *params) {
    double rijm,l;
    double *phi = (double *) malloc(6*sizeof(double));

    double *p = (double *) params;
    
    rijm = *d-1.0;
    phi[0] = rijm*rijm;
    
    if ( p[0] == 1 ) {
        l = p[1];
        phi[1] = (2/l)*(*d)*rijm;
    }
    
    return phi;
}






/*
int qadLJ(struct MCState *mcs,unsigned long int *nm, double *d) {
        
	unsigned long int ii,jj,dj,ind;
	double rij1,rij3,rij6,rij7,rij12,rij13;
        
	#pragma omp single
	{
	if ( nm == NULL) {
		mcs->nm = gsl_rng_uniform_int(mcs->rangen,mcs->N);
	}
	else {
		mcs->nm = *nm;
	}
	if ( d == NULL ) {
		mcs->rn = gsl_rng_uniform(mcs->rangen);
	}
	else {
		mcs->rn = *d;
	}

	mcs->md = (mcs->rn - 0.5)*2*mcs->maxStep;
	mcs->rTrial[mcs->nm] = mcs->r[mcs->nm] + mcs->md;
	}
	
	if (fabs(mcs->rTrial[mcs->nm]) > mcs->l/2.0) {
		#pragma omp single
		{
		mcs->dAcc[1]++;
		}
	}
	else {
		#pragma omp single
		{
		mcsdE6 = 0;
		mcsdE12 = 0;
		mcsdVir6 = 0;
		mcsdVir12 = 0;
		}
		
		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) reduction(+:mcsdE6,mcsdE12,mcsdVir6,mcsdVir12) nowait
		for (ii = 0; ii < mcs->nm; ii++) {
			dj = mcs->nm - ii - 1;
			ind = mcs->indind[ii][dj];
			mcs->rijTrial[ind] = mcs->rij[ind] + mcs->md;
			if ( mcs->nbn >= 0 && dj >= mcs->nbn ) {
				mcs->e6ijTrial[ind] = 0;
				mcs->e12ijTrial[ind] = 0;
				mcs->vir6ijTrial[ind] = 0;
				mcs->vir12ijTrial[ind] = 0;
			}
			else {
				rij1 = mcs->rijTrial[ind];
				rij3 = rij1*rij1*rij1;
				rij6 = 1/(rij3*rij3);
				rij7 = rij6/rij1;
				rij12 = rij6*rij6;
				rij13 = rij12/rij1;
				mcs->e6ijTrial[ind] = 4*rij6;
				mcs->e12ijTrial[ind] = 4*rij12;
				mcs->vir6ijTrial[ind] = 24/mcs->l*rij6;
				mcs->vir12ijTrial[ind] = 48/mcs->l*rij12;
			}
			mcsdE6  = mcsdE6  - mcs->e6ij[ind]  + mcs->e6ijTrial[ind];
			mcsdE12 = mcsdE12 - mcs->e12ij[ind] + mcs->e12ijTrial[ind];
			mcsdVir6 = mcsdVir6 - mcs->vir6ij[ind] + mcs->vir6ijTrial[ind];
			mcsdVir12 = mcsdVir12 - mcs->vir12ij[ind] + mcs->vir12ijTrial[ind];
		}

		#pragma omp for private(ii,jj,dj,ind,rij1,rij3,rij6,rij7,rij12,rij13) reduction(+:mcsdE6,mcsdE12,mcsdVir6,mcsdVir12)
		for (ii = mcs->nm + 1; ii < mcs->N; ii++) {
			dj = ii - mcs->nm - 1;
			ind = mcs->indind[mcs->nm][dj];
			mcs->rijTrial[ind] = mcs->rij[ind] - mcs->md;
			if ( mcs->nbn >= 0 && dj >= mcs->nbn ) {
				mcs->e6ijTrial[ind]    = 0;
				mcs->e12ijTrial[ind]   = 0;
				mcs->vir6ijTrial[ind]  = 0;
				mcs->vir12ijTrial[ind] = 0;
			}
			else {
				rij1 = mcs->rijTrial[ind];
				rij3 = rij1*rij1*rij1;
				rij6 = 1/(rij3*rij3);
				rij7 = rij6/rij1;
				rij12 = rij6*rij6;
				rij13 = rij12/rij1;
				mcs->e6ijTrial[ind] = 4*rij6;
				mcs->e12ijTrial[ind] = 4*rij12;
				mcs->vir6ijTrial[ind] = 24/mcs->l*rij6;
				mcs->vir12ijTrial[ind] = 48/mcs->l*rij12;
			}
			mcsdE6  = mcsdE6  - mcs->e6ij[ind]  + mcs->e6ijTrial[ind];
			mcsdE12 = mcsdE12 - mcs->e12ij[ind] + mcs->e12ijTrial[ind];
			mcsdVir6 = mcsdVir6 - mcs->vir6ij[ind] + mcs->vir6ijTrial[ind];
			mcsdVir12 = mcsdVir12 - mcs->vir12ij[ind] + mcs->vir12ijTrial[ind];
		}
		
		#pragma omp barrier  // THIS BARRIER IS CRITICAL TO CALCULATE dE ACCURATELY (I DO NOT UNDERSTAND WHY THAT SHOULD BE THE CASE)
		
		#pragma omp single
		{
		mcs->dE12 = mcsdE12;
		mcs->dE6 = mcsdE6;
		mcs->dE = mcs->dE12 - mcs->dE6;
		}
		
		if (mcs->dE > 0 ) {
		#pragma omp single
		{
			ran = gsl_rng_uniform(mcs->rangen);
			mcs->bf = exp(-mcs->dE/mcs->T);
		}
		}
		
		if (mcs->dE <= 0 || mcs->bf > ran) {
			#pragma omp single
			{
			mcs->dAcc[0]++;
			mcs->E      += mcs->dE;
			mcs->E6     += mcs->dE6;
			mcs->E12    += mcs->dE12;
			mcs->Vir6   += mcs->dVir6;
			mcs->Vir12  += mcs->dVir12;
			mcs->Vir    += mcs->dVir12 - mcs->dVir6;
			mcs->r[mcs->nm] = mcs->rTrial[mcs->nm];
			}
			
			#pragma omp for private(jj,dj,ind) nowait
			for (jj = 0; jj < mcs->nm; jj++) {
				dj = mcs->nm - jj - 1;
				ind = mcs->indind[jj][dj];
				mcs->rij[ind] = mcs->rijTrial[ind];
				mcs->e6ij[ind] = mcs->e6ijTrial[ind];
				mcs->e12ij[ind] = mcs->e12ijTrial[ind];
				mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
			}
			#pragma omp for private(jj,dj,ind)
			for (jj = mcs->nm + 1; jj < mcs->N; jj++) {
				dj = jj - mcs->nm - 1;
				ind = mcs->indind[mcs->nm][dj];
				mcs->rij[ind] = mcs->rijTrial[ind];
				mcs->e6ij[ind] = mcs->e6ijTrial[ind];
				mcs->e12ij[ind] = mcs->e12ijTrial[ind];
				mcs->vir6ij[ind] = mcs->vir6ijTrial[ind];
				mcs->vir12ij[ind] = mcs->vir12ijTrial[ind];
			}
			
                        // Short-cut calculation of density and g(x)
			qagrho(mcs);
			
		}
		else {
		#pragma omp single
		{
			mcs->dAcc[1]++;
		}
		}
	}
	
        // Accumulate density and g(x)
	ugrho(mcs);
	
        fflush(stdout);
	#pragma omp barrier

	return 0;
}



int qavHarmonic(struct MCState *mcs) {
	unsigned long int ii,ind;
	double dQ,ran,lRat1,lRat3,lRat6,lRat7,lRat12,lRat13;
	
	#pragma omp single
	{
	mcs->dl = (mcs->rn-0.5)*2*mcs->maxdl;
	lRat1 = (mcs->l+mcs->dl)/mcs->l;
	lRat3 = lRat1*lRat1*lRat1;
	lRat6 = 1/(lRat3*lRat3);
	lRat12 = lRat6*lRat6;
	mcs->E6Trial = lRat6*mcs->E6;
	mcs->E12Trial = lRat12*mcs->E12;
	mcs->dE = mcs->E12Trial - mcs->E6Trial - mcs->E;
	mcs->bf   = exp(-(mcs->dE + mcs->P*mcs->dl)/mcs->T + mcs->N*log(lRat1));

	if (mcs->bf < 1.0 ) {
		ran = gsl_rng_uniform(mcs->rangen);
	}
	}
	
	if (mcs->bf >= 1.0 || mcs->bf > ran) {
		#pragma omp single
		{
		mcs->vAcc[0]++;
		mcs->E6  = mcs->E6Trial;
		mcs->E12 = mcs->E12Trial;
		mcs->E      = mcs->E12 - mcs->E6;
		
		mcs->l = mcs->l + mcs->dl;
		
		lRat7 = lRat6/lRat1;
		lRat13 = lRat12/lRat1;
		mcs->Vir6 = lRat7*mcs->Vir6;
		mcs->Vir12 = lRat13*mcs->Vir12;
		mcs->Vir    = mcs->N*mcs->T/mcs->l + mcs->Vir12 - mcs->Vir6;
		}
		
		#pragma omp single private(ii) nowait
		for (ii = 0; ii < mcs->N; ii++) {
			mcs->r[ii] = lRat1*mcs->r[ii];
		}
		
		#pragma omp for private(ind)
		for (ind = 0; ind < mcs->numPairs; ind++) { 
			mcs->rij[ind]     = lRat1*mcs->rij[ind];
			mcs->e6ij[ind]    = lRat6*mcs->e6ij[ind];
			mcs->e12ij[ind]   = lRat12*mcs->e12ij[ind];
			mcs->vir6ij[ind]  = lRat7*mcs->vir6ij[ind];
			mcs->vir12ij[ind] = lRat13*mcs->vir12ij[ind];
		}
		
		fgrho(mcs);
	}
	
	else {
		#pragma omp single
		{
		mcs->vAcc[1]++;
		}
	}
	ugrho(mcs);
	fflush(stdout);

	return 0;
}
*/	

