#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



typedef struct PhaseState {
  double x, y,vx, vy,dx, dy,dvx , dvy,dax , day;
} PhaseState;

typedef struct PhaseStateSimple {
  double x, y,vx, vy;
} PhaseStateSimple;

typedef struct MEGNO_Auxilary_Variables {
  double W, Y, megno;
} MEGNO_Auxilary_Variables;

typedef struct Simulation {
	PhaseState * test_particle;
	PhaseStateSimple * inner_planet, * outer_planet;
	MEGNO_Auxilary_Variables * megno_aux;
	double mu1,mu2;
	double t;
} Simulation;

void simulation_step(Simulation * sim, double t, double dt);
void kepler_advance(Simulation * sim, double _dt);
void initialize_particle( PhaseState* particle, double n0, double lambda, double  ecc, double pomega);
void kepler_2D_advance_simple(  PhaseStateSimple* restrict particle ,double _dt);
void kepler_2D_advance(  PhaseState* restrict particle ,double _dt);
void intialize_simulation( Simulation * sim, double mu1, double mu2,
 double n1, double lambda1, double  ecc1, double pomega1,
 double n2, double lambda2, double  ecc2, double pomega2,
 double ntp, double lambdatp, double  ecctp, double pomegatp
 );
 void free_simulation( Simulation * simulation );
 
 double IntegrateSimulation(Simulation * sim, const double tFin, const double dt);
