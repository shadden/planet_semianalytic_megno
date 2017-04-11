#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)

typedef struct ResonanceData {
	int IncludeZeroth;
	int Nres,MaxOrder,MaxJ;
	int* ResonanceIndices;
	double* ResonanceCoefficients;
} ResonanceData;

typedef struct ActionAnglePhaseState {
	double L,l,Y,X;
	double dL,dl,dY,dX;
	double dLdot,dldot,dYdot,dXdot;
} ActionAnglePhaseState;

typedef struct SimulationParameters {
	double mu1,mu2,n1,n2,e1,e2,lambda2, varpi2;

	double alphaIn,alphaOut;
	double fSecIn,gSecIn;
	double fSecOut,gSecOut;

} SimulationParameters;

typedef struct PhaseState {
  double x, y,vx, vy,dx, dy,dvx , dvy,dax , day;
} PhaseState;

typedef struct MEGNO_Auxilary_Variables {
  double W, Y, megno;
} MEGNO_Auxilary_Variables;

typedef struct SecularSystem {
  double sx1,sy1,sx2,sy2;
  double * freqs;
  double * evrs ;
} SecularSystem;

typedef struct ActionAngleSimulation {
	ActionAnglePhaseState state;
	SimulationParameters parameters;
	MEGNO_Auxilary_Variables megno;
	ResonanceData rIn,rOut;
	double t;
	double dt;
} ActionAngleSimulation;


void initialize_pars(SimulationParameters* pars,double mu1,double mu2,double n1,double n2,double e1,double e2,double lambda2,double varpi2);
void intialize_megno_vars( MEGNO_Auxilary_Variables* megno);
void kepler_2D_advance(  PhaseState* restrict particle ,double _dt);
void two_circular_perturbers_advance(double const mu1, double const mu2, double const Omega2,  PhaseState* restrict particle , double t, double _dt);
void rotaion_advance(  PhaseState* restrict particle , double _dt);
void outer_circular_perturber_advance(double const mu1,  PhaseState* restrict particle ,double _dt);
void compute_variational_accs(PhaseState* restrict particle, double const mu1, double const mu2, double const Omega2, double t);
void update_megno_eqations(PhaseState* restrict particle , MEGNO_Auxilary_Variables* megno ,double t,double dt);
void intialize_particle( PhaseState* particle, double period, double  ecc);
void intialize_megno_vars( MEGNO_Auxilary_Variables* megno);

double secularF2(double alpha);
double secularF10(double alpha);
