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
#include <stdio.h>
#include <assert.h>
#include "shared.h"
#include "LaplaceCoefficients.h"

double laplace(double s, int i, int j, double a);
double FirstOrder(int k,int epower,double a);
double SecondOrder(int k,int epower,double a);
double FirstOrderI1O0(int k,double a);
double FirstOrderI0O1(int k,double a);
double SecondOrderI2O0(int k,double a);
double SecondOrderI1O1(int k,double a);
double SecondOrderI0O2(int k,double a);
double GeneralOrderCoefficient(int res_j, int order, int epower,double a);

int binomialCoeff(int n, int k)
{
  // Base Cases
  if (k==0 || k==n)
    return 1;
 
  // Recur
  return  binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k);
}
int factorial(int f)
{
 	assert(f>=0);
    if ( f == 0 ) 
        return 1;
    return(f * factorial(f - 1));
}


void initialize_ResonanceData( ResonanceData * r,int inclue0th){
	r->Nres =0;
	r->MaxOrder=0;
	r->MaxJ=0;
	r->ResonanceIndices=NULL;
	r->ResonanceCoefficients=NULL;
	r-> IncludeZeroth = inclue0th;
}

void free_ResonanceData( ResonanceData * r){
	free(r->ResonanceIndices);
	free(r->ResonanceCoefficients);
}

void AddResonance(ResonanceData* r, int res_j, int order, int epower,double a){
	assert(order >= epower && epower >=0);

	double coeff;

	// Only resonances up to MAX_ORDER allowed
	assert(order<= MAX_ORDER);
	assert(res_j < MAX_J);

	coeff=GeneralOrderCoefficient(res_j,order,epower,a);

	r->Nres +=1;

	if (r->MaxOrder < order){
		r->MaxOrder = order;
	}
	if (r->MaxJ < res_j){
		r->MaxJ = res_j;
	}

	r->ResonanceIndices = (int *) realloc(r->ResonanceIndices,  (r->Nres) * 3 * sizeof(int));
	r->ResonanceCoefficients = (double *) realloc(r->ResonanceCoefficients,  (r->Nres) * ( MAX_ORDER + 1 ) *  sizeof(double));
	
	const int NLast = r->Nres-1;

	*(r->ResonanceIndices + 3*NLast) = res_j;
	*(r->ResonanceIndices + 3*NLast + 1) = order;
	*(r->ResonanceIndices + 3*NLast + 2) = epower;
	
	*( r->ResonanceCoefficients + ( MAX_ORDER + 1 )*NLast + epower ) = coeff;
	
	
}

double secularF2(double alpha){
	assert(alpha<1);
	return 0.25 * laplace(0.5,0,1,alpha) + 0.125 * laplace(0.5,0,2,alpha);
};
double secularF10(double alpha){
		assert(alpha<1);
		return 0.5 * laplace(0.5,1,0,alpha) - 0.5 * laplace(0.5,1,1,alpha) - 0.25 * laplace(0.5,1,2,alpha);
};

// Page 247 of M&D '99

// Leading order component of M & D Equation 6.38

double NCOd0(int a, int b, int c){
	if(c == 0) return 1.0;
	if(c==1) return b - 0.5 * a;
	// Else
	double cc = (double) c;
	double nc1,nc2;
	nc1 = NCOd0(a, b+1,c-1);
	nc2 = NCOd0(a, b+2,c-2);
	return 0.25 * (2 * (2*b - a) * nc1 + (b - a) * nc2 ) / cc;
}
double NewcombOperator(int a, int b, int c,int d){
	assert(c==0 || d==0);
	if (d==0){
		return NCOd0(a,b,c);
	}
	else{
		return NCOd0(a,-b,d);
	}
}
double HansenCoefficient(int a, int b, int c){
	int alpha =  c-b > 0 ? c-b:0;
	int beta =  b-c > 0 ? b-c:0;
	return NewcombOperator(a,b,alpha,beta);
}

double GeneralOrderCoefficient(int res_j, int order, int epower,double a){

	int j[7];
	j[0]=0; // ignore to match M&D indexing
	j[1] = res_j;
	j[2] = order - res_j ;
	j[3] = epower - order;
	j[4] = -1 * epower;
	j[5] = 0;
	j[6] = 0; 
	int q = j[4];
	int q1 = -1 * j[3];
	int Nmax = order;
	
	double coeff =0;
	for (int l=0; l<= Nmax; l++){
	int sgn = l%2 ? -1:1;
	double fact = (double) factorial(l);
	double sum = 0;
	for (int k=0; k<= l; k++){
		double ncIn;
		double ncOut;
		int sgn2 = k%2 ? -1:1;
		int binom = sgn2 * binomialCoeff(l, k);
		int jj = j[2] + q;
		if (jj<0) jj=-1*jj;
		ncIn = HansenCoefficient(k,-j[2]-j[4],-j[2]);
		ncOut = HansenCoefficient(-1-k,j[1]+j[3],j[1]);
		sum += binom * ncIn*ncOut* laplace(0.5,jj,l,a);
	}
	coeff += sgn * sum / fact;
	}
	return coeff;

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




//	**** DEPRICATED **** //

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

