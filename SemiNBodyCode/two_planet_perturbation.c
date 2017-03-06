#include "SemiNBody.h"


static inline double random_real(){
        return -1+2.*((float)rand())/RAND_MAX;
}


void interaction_advance(Simulation * sim, double _dt){
	const PhaseState p1 = *(sim->test_particle);
	const double mu1 = sim->mu1;
	const double mu2 = sim->mu2;	

	// planet-1 separation
	PhaseStateSimple * inner_planet = sim->inner_planet;
	const double x1 = inner_planet->x;
	const double y1 = inner_planet->y;
	const double dx1 = x1 - p1.x ;
	const double dy1 = y1 - p1.y ;
	const double rho1sq = dx1 * dx1  + dy1 * dy1 ;
	const double rho1 = sqrt( rho1sq );
	const double rho1_inv3 = 1./(rho1sq * rho1);
	const double rho1_inv5 = rho1_inv3/rho1sq;
	
	// planet-2 separation

	PhaseStateSimple * outer_planet = sim->outer_planet;
	const double x2 = outer_planet->x;
	const double y2 = outer_planet->y;
	const double dx2 = x2 - p1.x ;
	const double dy2 = y2 - p1.y ;
	const double rho2sq = dx2 * dx2  + dy2 * dy2 ;
	
	const double rho2 = sqrt( rho2sq );
	
	const double rho2_inv3 = 1./(rho2sq * rho2);
	const double rho2_inv5 = rho2_inv3/rho2sq;
	
	const double r2sq = x2*x2 +y2*y2;
	const double r2cubed_inv = 1. / pow(r2sq,1.5);
	
	// test-particle r	
	const double rsq = p1.x * p1.x + p1.y * p1.y;
	const double r = sqrt(rsq);
	const double r_inv3 = 1. / rsq / r; 
	const double r_inv5 = r_inv3 / rsq;
	const double r_inv7 = r_inv5 / rsq;
	const double r_dot_r1 = p1.x*x1 + p1.y*y1;

	const double F1x = dx1 * rho1_inv3 + 3 * r_dot_r1 * p1.x * r_inv5 - x1 * r_inv3;
	const double F1y = dy1 * rho1_inv3 + 3 * r_dot_r1 * p1.y * r_inv5 - y1 * r_inv3;
	
	const double F2x = dx2 * rho2_inv3 -  x2  * r2cubed_inv;
	const double F2y = dy2 * rho2_inv3 -  y2  * r2cubed_inv;
	
	PhaseState * particle = sim->test_particle;
	(*particle).vx += _dt * ( mu1 * F1x + mu2 * F2x );
	(*particle).vy += _dt * ( mu1 * F1y + mu2 * F2y );
	

	// Variationals
	
	const double DF1xx =  3  * dx1*dx1 * rho1_inv5 - rho1_inv3 \
						- 15 * p1.x * p1.x * r_dot_r1 * r_inv7\
					 	+ 6  * p1.x * x1 * r_inv5\
					 	+ 3  * r_dot_r1 * r_inv5;

	const double DF1yy =  3  * dy1*dy1 * rho1_inv5 - rho1_inv3 \
						- 15 * p1.y * p1.y * r_dot_r1 * r_inv7\
					 	+ 6  * p1.y * y1 * r_inv5\
					 	+ 3  * r_dot_r1 * r_inv5;

	const double DF1xy =  3  * dx1*dy1 * rho1_inv5 \
						 -15 * p1.x * p1.y * r_dot_r1 * r_inv7\
						 +3  * (x1*p1.y+y1*p1.x) * r_inv5;

	const double DF2xx = 3*dx2*dx2*rho2_inv5-rho2_inv3;
	const double DF2xy = 3*dx2*dy2*rho2_inv5;
	const double DF2yy = 3*dy2*dy2*rho2_inv5-rho2_inv3;

	// inner planet variational 
	const double Dax = mu1 * (DF1xx * p1.dx + DF1xy * p1.dy ) + mu2 * (DF2xx * p1.dx + DF2xy * p1.dy );
	const double Day = mu1 * (DF1xy * p1.dx + DF1yy * p1.dy ) + mu2 * (DF2xy * p1.dx + DF2yy * p1.dy );

	(*particle).dvx += _dt * Dax;
	(*particle).dvy += _dt * Day;


}

void kepler_advance(Simulation * sim, double _dt){
	kepler_2D_advance(sim->test_particle,_dt);
	kepler_2D_advance_simple(sim->inner_planet,_dt);
	kepler_2D_advance_simple(sim->outer_planet,_dt);
}

void initialize_megno_vars( MEGNO_Auxilary_Variables* megno){

        (*megno).Y=0;
        (*megno).W=0;
        (*megno).megno=0;

}

void initialize_particle( PhaseState* particle, double n0, double lambda, double  ecc, double pomega){

	const double a0 = pow(n0*n0, -1/3. );
	
	// Positions and velocities in orbit reference frame	
	const double x0 = a0 * (1 - ecc);
	const double y0 = 0;
	const double vx0 = 0;
	const double vy0 = n0*a0*sqrt(1+ecc)/sqrt(1-ecc) ;

	// rotate by periapse angle
	const double sw = sin(pomega);
	const double cw = cos(pomega);
	(*particle).x = cw * x0 - sw * y0;
	(*particle).y = sw * x0 + cw * y0;
	(*particle).vx= cw * vx0 - sw * vy0;
	(*particle).vy= sw * vx0 + cw * vy0;
	
	// lambda = M + pomega
	const double M0 = lambda - pomega;
	const double dt0 =	M0 / n0;
	kepler_2D_advance(particle, dt0);
	
	double dvec[4];
	double normsq=0;
	
	// generate random reals	
	for(int i=0; i < 4; i++){
		dvec[i] = random_real();
		normsq += dvec[i]*dvec[i];
	}
	
	(*particle).dx = 	1.e-15*dvec[0]/sqrt(normsq);
	(*particle).dy = 	1.e-15*dvec[1]/sqrt(normsq);
	(*particle).dvx = 1.e-15*dvec[2]/sqrt(normsq);
	(*particle).dvy = 1.e-15*dvec[3]/sqrt(normsq);

}

void initialize_particle_simple( PhaseStateSimple* particle, double n0, double lambda, double  ecc, double pomega){

	const double a0 = pow(n0*n0, -1/3. );
	
	// Positions and velocities in orbit reference frame	
	const double x0 = a0 * (1 - ecc);
	const double y0 = 0;
	const double vx0 = 0;
	const double vy0 = n0*a0*sqrt(1+ecc)/sqrt(1-ecc) ;

	// rotate by periapse angle
	const double sw = sin(pomega);
	const double cw = cos(pomega);
	(*particle).x = cw * x0 - sw * y0;
	(*particle).y = sw * x0 + cw * y0;
	(*particle).vx= cw * vx0 - sw * vy0;
	(*particle).vy= sw * vx0 + cw * vy0;
	
	// lambda = M + pomega
	const double M0 = lambda - pomega;
	const double dt0 =	M0 / n0;
	kepler_2D_advance_simple(particle, dt0);
}

void intialize_simulation( Simulation * sim, double mu1, double mu2,
 double n1, double lambda1, double  ecc1, double pomega1,
 double n2, double lambda2, double  ecc2, double pomega2,
 double ntp, double lambdatp, double  ecctp, double pomegatp
 ){
	sim->mu1=mu1;
	sim->mu2=mu2;
	sim->t=0.0;

	// allocate memory
	sim->test_particle = (PhaseState *) malloc(sizeof(PhaseState));
	sim->inner_planet = (PhaseStateSimple *) malloc(sizeof(PhaseStateSimple));
	sim->outer_planet = (PhaseStateSimple *) malloc(sizeof(PhaseStateSimple));
	sim->megno_aux = (MEGNO_Auxilary_Variables *) malloc(sizeof(MEGNO_Auxilary_Variables));
 
 	// initialize particles
 	initialize_particle(sim->test_particle,ntp,lambdatp,ecctp,pomegatp);
 	initialize_particle_simple(sim->inner_planet,n1,lambda1,ecc1,pomega1);
 	initialize_particle_simple(sim->outer_planet,n2,lambda2,ecc2,pomega2);
 	initialize_megno_vars(sim->megno_aux);
 
}

void free_simulation( Simulation * simulation ){
 free(simulation->test_particle);
 free(simulation->inner_planet);
 free(simulation->outer_planet);
 free(simulation->megno_aux); 
}

void update_megno_equations(Simulation * restrict sim , double t, double dt){
	const PhaseState p1 = *(sim->test_particle);
	MEGNO_Auxilary_Variables m1 = *(sim->megno_aux);
	double deltaSq = p1.dx*p1.dx + p1.dy*p1.dy + p1.dvx*p1.dvx + p1.dvy*p1.dvy;
	double deltad_delta = p1.dx*p1.dvx + p1.dy*p1.dvy + p1.dax*p1.dvx + p1.day*p1.dvy;

	MEGNO_Auxilary_Variables * megno = sim->megno_aux;
	(*megno).Y = m1.Y +  deltad_delta / deltaSq * t * dt;
	(*megno).W = m1.W + dt * 2 * m1.Y / t ;
	(*megno).megno = m1.W/t;
}

void compute_variational_accs(Simulation * sim){

	// test-particle planet-1 separation
	PhaseState * particle=sim->test_particle;
	const PhaseState p1 = *(sim->test_particle);
	const double mu1 = sim->mu1;
	const double mu2 = sim->mu2;	

	// planet-1 separation
	PhaseStateSimple * inner_planet = sim->inner_planet;
	const double x1 = inner_planet->x;
	const double y1 = inner_planet->y;
	const double dx1 = x1 - p1.x ;
	const double dy1 = y1 - p1.y ;
	const double rho1sq = dx1 * dx1  + dy1 * dy1 ;
	const double rho1 = sqrt( rho1sq );
	const double rho1_inv3 = 1./(rho1sq * rho1);
	const double rho1_inv5 = rho1_inv3/rho1sq;
	
	// planet-2 separation

	PhaseStateSimple * outer_planet = sim->outer_planet;
	const double x2 = outer_planet->x;
	const double y2 = outer_planet->y;
	const double dx2 = x2 - p1.x ;
	const double dy2 = y2 - p1.y ;
	const double rho2sq = dx2 * dx2  + dy2 * dy2 ;
	
	const double rho2 = sqrt( rho2sq );
	
	const double rho2_inv3 = 1./(rho2sq * rho2);
	const double rho2_inv5 = rho2_inv3/rho2sq;
	
	const double r2sq = x2*x2 +y2*y2;
	const double r2cubed_inv = 1. / pow(r2sq,1.5);
	
	// test-particle r	
	const double rsq = p1.x * p1.x + p1.y * p1.y;
	const double r = sqrt(rsq);
	const double r_inv3 = 1. / rsq / r; 
	const double r_inv5 = r_inv3 / rsq;
	const double r_inv7 = r_inv5 / rsq;
	const double r_dot_r1 = p1.x*x1 + p1.y*y1;
	
	const double DFCENxx = 3 * p1.x * p1.x * r_inv5 - r_inv3;
	const double DFCENxy = 3 * p1.x * p1.y * r_inv5 ;
	const double DFCENyy = 3 * p1.y * p1.y * r_inv5 - r_inv3 ;
	// Variationals
	const double DF1xx =  3  * dx1*dx1 * rho1_inv5 - rho1_inv3 \
						- 15 * p1.x * p1.x * r_dot_r1 * r_inv7\
					 	+ 6  * p1.x * x1 * r_inv5\
					 	+ 3  * r_dot_r1 * r_inv5;

	const double DF1yy =  3  * dy1*dy1 * rho1_inv5 - rho1_inv3 \
						- 15 * p1.y * p1.y * r_dot_r1 * r_inv7\
					 	+ 6  * p1.y * y1 * r_inv5\
					 	+ 3  * r_dot_r1 * r_inv5;

	const double DF1xy =  3  * dx1*dy1 * rho1_inv5 \
						 -15 * p1.x * p1.y * r_dot_r1 * r_inv7\
						 +3  * (x1*p1.y+y1*p1.x) * r_inv5;

	const double DF2xx = 3*dx2*dx2*rho2_inv5-rho2_inv3;
	const double DF2xy = 3*dx2*dy2*rho2_inv5;
	const double DF2yy = 3*dy2*dy2*rho2_inv5-rho2_inv3;

	
	(*particle).dax =   (DFCENxx + mu1*DF1xx+ mu2*DF2xx) * p1.dx + (DFCENxy + mu1*DF1xy+ mu2*DF2xy) * p1.dy ;
	(*particle).day =   (DFCENxy + mu1*DF1xy+ mu2*DF2xy) * p1.dx + (DFCENyy + mu1*DF1yy+ mu2*DF2yy) * p1.dy ;
	
}

void simulation_step(Simulation * sim, double t, double dt){
	// update positions, velocities
	interaction_advance(sim,0.5*dt);
	kepler_advance(sim,dt);
	interaction_advance(sim,0.5*dt);
	
	//update variationals, megno
	compute_variational_accs(sim);
	update_megno_equations(sim,t+dt,dt);
	
	
}

double IntegrateSimulation(Simulation * sim, const double tFin, const double dt){
	const int Nstep = (int) ceil( tFin / dt);

	for(int i=0; i<Nstep; i++){
		simulation_step(sim,i*dt,dt);
	}
	return sim->megno_aux->megno;
}


void IntegrateSimulationToTime(Simulation * sim, const double tFin, const double dt){
	const double current_t = sim->t;
	const int Nstep = (int) ceil( (tFin-current_t) / dt);

	for(int i=0; i<Nstep; i++){
		simulation_step(sim,i*dt,dt);
	}
	sim->t = current_t +  Nstep * dt;
}
