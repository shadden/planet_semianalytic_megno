#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)


typedef struct PhaseState {
  double x, y,vx, vy,dx, dy,dvx , dvy,dax , day;
} PhaseState;

typedef struct MEGNO_Auxilary_Variables {
  double W, Y, megno;
} MEGNO_Auxilary_Variables;

void intialize_megno_vars( MEGNO_Auxilary_Variables* megno);
void kepler_2D_advance(  PhaseState* restrict particle ,double _dt);
void two_circular_perturbers_advance(double const mu1, double const mu2, double const Omega2,  PhaseState* restrict particle , double t, double _dt);
void rotaion_advance(  PhaseState* restrict particle , double _dt);
void outer_circular_perturber_advance(double const mu1,  PhaseState* restrict particle ,double _dt);
void compute_variational_accs(PhaseState* restrict particle, double const mu1, double const mu2, double const Omega2, double t);
void update_megno_eqations(PhaseState* restrict particle , MEGNO_Auxilary_Variables* megno ,double t,double dt);
void intialize_particle( PhaseState* particle, double period, double  ecc);
void intialize_megno_vars( MEGNO_Auxilary_Variables* megno);
