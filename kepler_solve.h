#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)



// typedef struct {
//   double x, y, vx, vy, dx, dy ,dvx , dvy;
// } PhaseState;

typedef struct PhaseState {
  double x, y, vx, vy, dx, dy ,dvx , dvy;
} PhaseState;

void kepler_2D_advance(  PhaseState* restrict particle ,double _dt);
