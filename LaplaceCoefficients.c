/*Many thanks for Jack Wisdom at MIT for sharing the code "laplace()" below with us. */
/* the longitude of transit = 0 */
/* Please cite Agol & Deck (2015) if you make use of this code in your research.*/
#define LAPLACE_EPS 1.0e-10

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define TWOPI 6.283185307179586476925287

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "Resonances.h"

double laplace(double s, int i, int j, double a);
double FirstOrder(int k,int epower,double a);
double SecondOrder(int k,int epower,double a);
double FirstOrderI1O0(int k,double a);
double FirstOrderI0O1(int k,double a);
double SecondOrderI2O0(int k,double a);
double SecondOrderI1O1(int k,double a);
double SecondOrderI0O2(int k,double a);


void initialize_ResonanceData( ResonanceData * r){
	r->Nres =0;
	r->MaxOrder=0;
	r->ResonanceIndices=NULL;
	r->ResonanceCoefficients=NULL;
}

void free_ResonanceData( ResonanceData * r){
	free(r->ResonanceIndices);
	free(r->ResonanceCoefficients);
}

void AddResonance(ResonanceData* r, int res_j, int order, int epower,double a){
	assert(order >= epower && epower >=0);

	double coeff;

	// Only resonances up to 2nd order supported currently
	assert(order<= MAX_ORDER);

	if(order==1){
		coeff = FirstOrder(res_j,epower,a);
	}else{
		coeff = SecondOrder(res_j,epower,a);
	}

	r->Nres +=1;

	if (r->MaxOrder < order){
		r->MaxOrder = order;
	}

	r->ResonanceIndices = (int *) realloc(r->ResonanceIndices,  (r->Nres) * 3 * sizeof(int));
	r->ResonanceCoefficients = (double *) realloc(r->ResonanceCoefficients,  (r->Nres) * ( MAX_ORDER + 1 ) *  sizeof(double));
	
	const int NLast = r->Nres-1;

	*(r->ResonanceIndices + 3*NLast) = res_j;
	*(r->ResonanceIndices + 3*NLast + 1) = order;
	*(r->ResonanceIndices + 3*NLast + 2) = epower;
	
	*( r->ResonanceCoefficients + ( MAX_ORDER + 1 )*NLast + epower ) = coeff;
	
	
}

double FirstOrder(int k,int epower,double a){
	if(epower==1){
		return FirstOrderI1O0(k,a);
	} else {
		return FirstOrderI0O1(k,a);
	}
}

double SecondOrder(int k,int epower,double a){
	switch(epower){
		case 2:
			return SecondOrderI2O0(k,a);
		case 1:
			return SecondOrderI1O1(k,a);
		default:
			return SecondOrderI0O2(k,a);
	}
}

double FirstOrderI1O0(int k,double a) {
	return -(k*laplace(0.5,k,0,a)) - laplace(0.5,k,1,a)/2.;
}

double FirstOrderI0O1(int k,double a) {
	return (-0.5 + k)*laplace(0.5,-1 + k,0,a) + laplace(0.5,-1 + k,1,a)/2.;
}

double SecondOrderI2O0(int k,double a) {
	return ((-5*k)/8. + k*k/2.)*laplace(0.5,k,0,a) + (-0.25 + k/2.)*laplace(0.5,k,1,a) + laplace(0.5,k,2,a)/8.;
}
double SecondOrderI1O1(int k,double a) {
	return (-0.5 + (3*k)/2. - k*k)*laplace(0.5,-1 + k,0,a) + (0.5 - k)*laplace(0.5,-1 + k,1,a) - laplace(0.5,-1 + k,2,a)/4.;
}
double SecondOrderI0O2(int k,double a) {
	return (0.25 - (7*k)/8. + k*k/2.)*laplace(0.5,-2 + k,0,a) + (-0.25 + k/2.)*laplace(0.5,-2 + k,1,a) + laplace(0.5,-2 + k,2,a)/8.;
}


/* Code due to Jack Wisdom */
/* compute Laplace coefficients and Leverrier derivatives
          j
     j   d     i
    a   ---   b (a)
          j    s
        da
   by series summation */



double laplace(double s, int i, int j, double a)
{
  double as, term, sum, factor1, factor2, factor3, factor4;
  int k,q, q0;

  as = a*a;

  if(i<0) i = -i;

  if(j<=i)     /* compute first term in sum */
    {
      factor4 = 1.0;
      for(k=0; k<j; k++)
        factor4 *= (i - k);
      sum = factor4;
      q0=0;
    }
  else
    {
       q0 = (j + 1 - i) / 2;
      sum = 0.0;
      factor4 = 1.0;
    }

  /* compute factors for terms in sum */

  factor1 = s;
  factor2 = s + i;
  factor3 = i + 1.0;
  for(q=1;q<q0;q++)   /* no contribution for q = 0 */
    {
      factor1 *= s + q;
      factor2 *= s + i + q;
      factor3 *= i + 1.0 + q;
    }

  term = as * factor1 * factor2 / (factor3 * q);

  /* sum series */

  while(term*factor4 > LAPLACE_EPS)
    {
      factor4 = 1.0;
      for(k=0;k<j;k++)
        factor4 *= (2*q + i - k);
      sum += term * factor4;
      factor1 += 1.0;
      factor2 += 1.0;
      factor3 += 1.0;
      q++;
      term *= as * factor1 * factor2 / (factor3 * q);

    }

  /* fix coefficient */

  for(k=0;k<i;k++)
    sum *= (s + ((double) k))/(((double) k)+1.0);

  if(q0 <= 0)
    sum *= 2.0 * pow(a, ((double) i));
  else
    sum *= 2.0 * pow(a, ((double) 2*q0 + i - 2));

  return(sum);
}

