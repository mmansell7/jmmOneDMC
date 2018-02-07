#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>


double * phiLJ(double *d, void *params) {
	double rij1,rij3,rij6,rij7,rij12,rij13,phi6,phi12,vir6,vir12;
	double *phi = (double *) malloc(4*sizeof(double));
	//printf("Check inside ####### d = %lf\n",d);

	rij1 = *d;
	rij3 = rij1*rij1*rij1;
	rij6 = 1/(rij3*rij3);
	rij7 = rij6/rij1;
	rij12 = rij6*rij6;
	rij13 = rij12/rij1;

	phi6 = 4*rij6;
	phi12 = 4*rij12;
	//vir6 = 24/l*rij6;
	//vir12 = 48/l*rij12;
	phi[0] = phi6;
	phi[1] = phi12;
	//phi[2] = vir6;
	//phi[3] = vir12;

	return phi;
}

