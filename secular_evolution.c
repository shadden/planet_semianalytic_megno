#include <stdlib.h>
#include <math.h>
#include "shared.h"
#include <stdio.h>

void get_secular_frequencies_and_vectors(double * freqs, double * eigen_vector_ratios, double mu1, double mu2, double fcoeff, double gcoeff, double n2, double alpha){

	double sqrt_alpha_inv = 1./sqrt(alpha);
	
	double s11 = mu2 * n2 * 2 * fcoeff * sqrt_alpha_inv;
	double s12 = mu2 * n2 * gcoeff * sqrt_alpha_inv;
	double s21 = mu1 * n2 * gcoeff ;
	double s22 = mu1 * n2 * 2 * fcoeff ;
	
	double trace = s11 + s22;
	double det = s11 * s22 - s12 * s21;
	// x^2 + trace * x + det = 0
	double disc_sq = trace * trace - 4 * det;

	*(freqs+0) = 0.5 *  (trace + sqrt(disc_sq));
	*(freqs+1) = 0.5 *  (trace - sqrt(disc_sq));
	
	*(eigen_vector_ratios+0) = -1 * s12 / ( s11 - *(freqs+0));
	*(eigen_vector_ratios+1) = -1 * s12 / ( s11 - *(freqs+1));
}

void initialize_SecularSystem(SecularSystem * s,SimulationParameters * pars){
	const double n1 = pars->n1;
	const double n2 = pars->n2;
	const double alpha = pow(n2/n1,2./3.);
	const double fcoeff = secularF2(alpha);
	const double gcoeff = secularF10(alpha);
	
	s->freqs = (double *) malloc(2 * sizeof(double));
	s->evrs  = (double *) malloc(2 * sizeof(double));
	get_secular_frequencies_and_vectors(s->freqs, s->evrs, pars->mu1, pars->mu1,fcoeff,gcoeff, n2, alpha);
	double x1=pars->e1;
	double y1=0;
	double x2=(pars->e2) * cos(pars->varpi2);
	double y2=-1*(pars->e2) * sin(pars->varpi2);
	
	double r1 = *(s->evrs);
	double r2 = *(s->evrs+1);
	double r_diff_inv = 1./(r1-r2);

	s->sx1 =  r_diff_inv * x1 - r_diff_inv * r2 * x2;
	s->sy1 =  r_diff_inv * y1 - r_diff_inv * r2 * y2;
	s->sx2 =  r_diff_inv * r1 * x2 - r_diff_inv * x1;
	s->sy2 =  r_diff_inv * r1 * y2 - r_diff_inv * y1;
}

void get_SecularSolution(double * ecc_components,SecularSystem * s,double t){
	
	double r1 = *(s->evrs);
	double r2 = *(s->evrs+1);
	double f1 = *(s->freqs);
	double f2 = *(s->freqs+1);

	double c1 = cos(f1 * t);
	double s1 = sin(f1 * t);
	double c2 = cos(f2 * t);
	double s2 = sin(f2 * t);
	
	// x1
	*(ecc_components) 	= r1 * ( (s->sx1) * c1 + (s->sy1) * s1) + r2 * ( (s->sx2) * c2 + (s->sy2) * s2);
	// y1
	*(ecc_components+1) = r1 * ( (s->sx1) * s1 - (s->sy1) * c1) + r2 * ( (s->sx2) * s2 - (s->sy2) * c2);

	// x2
	*(ecc_components+2) = ( (s->sx1) * c1 + (s->sy1) * s1) + ( (s->sx2) * c2 + (s->sy2) * s2);
	// y2
	*(ecc_components+3) = ( (s->sx1) * s1 - (s->sy1) * c1) + ( (s->sx2) * s2 - (s->sy2) * c2);
	
}


void free_SecularSystem(SecularSystem * s){
	free(s->freqs);
	free(s->evrs);
}