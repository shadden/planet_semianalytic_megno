#include "shared.h"

static inline double random_real(){
        return -1+2.*((float)rand())/RAND_MAX;
}
double MEGNO_Integration(double tfin, double dt, double period, double ecc, double mu1,double mu2, double Omega2){

	int Nstep = (int)(tfin/dt);

	PhaseState planet;
	intialize_particle(&planet,period,ecc);
	MEGNO_Auxilary_Variables megno;
	intialize_megno_vars( &megno);


	for(int i=0; i < Nstep; i++){
	// Both planets 
		two_circular_perturbers_advance( mu1 , mu2 , Omega2 , &planet , i*dt , dt*0.5 );	
		kepler_2D_advance(&planet, dt);
		rotaion_advance(&planet,dt);
		two_circular_perturbers_advance( mu1 , mu2 , Omega2 , &planet , i*dt + dt , dt*0.5 );

		compute_variational_accs(&planet,mu1,mu2,Omega2,i*dt + dt);
		update_megno_eqations( &planet , &megno ,i*dt + dt, dt);
	}

	return megno.megno;
	
}

void intialize_particle( PhaseState* particle, double period, double  ecc){

	const double a0 = pow(period*period, 1/3. );
	(*particle).x = a0 * (1 - ecc);
	(*particle).y = 0;
	(*particle).vx= 0;
	(*particle).vy= sqrt(1+ecc)/sqrt(a0) ;
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

void two_circular_perturbers_advance(double const mu1, double const mu2, double const Omega2,  PhaseState* restrict particle , double t, double _dt){
	const PhaseState p1 = *particle;
	
	// test-particle planet-1 separation
	const double dx1 = 1. - p1.x ;
	const double rho1sq = dx1 * dx1  + p1.y * p1.y ;
	const double rho1 = sqrt( rho1sq );
	const double rho1_inv3 = 1./(rho1sq * rho1);
	const double rho1_inv5 = rho1_inv3/rho1sq;
	// planet-2 location
	const double r2 = pow(Omega2, -2./3.);
	const double theta2 = Omega2 * t - t;
	const double x2 =  r2 * cos(theta2) ;
	const double y2 =  r2 * sin(theta2) ;	
	const double dx2 =  x2-p1.x ;
	const double dy2 =  y2-p1.y ;	
	const double rho2sq = dx2*dx2 + dy2*dy2;
	const double rho2 = sqrt(rho2sq);
	const double rho2_inv3 = 1./rho2sq / rho2;
	const double rho2_inv5 = rho2_inv3/rho2sq;
	// test-particle r	
	const double rsq = p1.x * p1.x + p1.y * p1.y;
	const double r = sqrt(rsq);
	const double r_inv3 = 1. / rsq / r; 
	const double r_inv5 = r_inv3 / rsq; 

	const double F1x = dx1 * rho1_inv3 + 3. * p1.x * p1.x * r_inv5 - r_inv3;
	const double F1y = -p1.y * rho1_inv3 + 3. * p1.x * p1.y * r_inv5;
	const double F2x = dx2 * rho2_inv3 -  x2  * Omega2 * Omega2;
	const double F2y = dy2 * rho2_inv3 -  y2  * Omega2 * Omega2;
	
	
	(*particle).vx += _dt * ( mu1 * F1x + mu2 * F2x );
	(*particle).vy += _dt * ( mu1 * F1y + mu2 * F2y );
	
	// Variationals
	const double DF1xx = 3.* dx1*dx1 * rho1_inv5 - rho1_inv3 \
					 	- 15.*p1.x*p1.x*p1.x * r_inv5 / rsq \
					 	+ 9*p1.x * r_inv5;
	const double DF1xy = -3. * dx1 * p1.y * rho1_inv5 \
						- 15.*p1.x * p1.x * p1.y * r_inv5 / rsq \
						+ 3*p1.y * r_inv5;
	const double DF1yy = ( 2*p1.y*p1.y - dx1*dx1 ) * rho1_inv5 + 3*p1.x*(p1.x*p1.x-4*p1.y*p1.y) * r_inv5 / rsq;
	const double DF2xx = 3*dx2*dx2*rho2_inv5-rho2_inv3;
	const double DF2xy = 3*dx2*dy2*rho2_inv5;
	const double DF2yy = 3*dy2*dy2*rho2_inv5-rho2_inv3;


	// inner planet variational 
	const double Dax = mu1 * (DF1xx * p1.dx + DF1xy * p1.dy ) + mu2 * (DF2xx * p1.dx + DF2xy * p1.dy );
	const double Day = mu1 * (DF1xy * p1.dx + DF1yy * p1.dy ) + mu2 * (DF2xy * p1.dx + DF2yy * p1.dy );

	(*particle).dvx += _dt * Dax;
	(*particle).dvy += _dt * Day;
	

}


void rotaion_advance(  PhaseState* restrict particle , double _dt){
	const PhaseState p1 = *particle;
	// eqns
	(*particle).x = cos(_dt)*p1.x + sin(_dt)*p1.y;
	(*particle).y = cos(_dt)*p1.y - sin(_dt)*p1.x;
	(*particle).vx = cos(_dt)*p1.vx + sin(_dt)*p1.vy;
	(*particle).vy = cos(_dt)*p1.vy - sin(_dt)*p1.vx;	
	// variationals
	(*particle).dx = cos(_dt)*p1.dx + sin(_dt)*p1.dy;
	(*particle).dy = cos(_dt)*p1.dy - sin(_dt)*p1.dx;
	(*particle).dvx = cos(_dt)*p1.dvx + sin(_dt)*p1.dvy;
	(*particle).dvy = cos(_dt)*p1.dvy - sin(_dt)*p1.dvx;	


}
void outer_circular_perturber_advance(double const mu1,  PhaseState* restrict particle ,double _dt){
	const PhaseState p1 = *particle;
	
	// test-particle planet-1 separation

	const double dx1 = 1. - p1.x ;
	const double rho1sq = dx1 * dx1  + p1.y * p1.y ;
	const double rho1 = sqrt( rho1sq );
	
	const double rsq = p1.x * p1.x + p1.y * p1.y;
	const double r = sqrt(rsq);
	
	const double Fx = dx1 / ( rho1sq * rho1) - 1.0;
	const double Fy = -p1.y /( rho1sq * rho1);

	(*particle).vx += _dt * ( mu1 * Fx  );
	(*particle).vy += _dt * ( mu1 * Fy  );

}
void update_megno_eqations(PhaseState* restrict particle , MEGNO_Auxilary_Variables* megno ,double t,double dt){
	const PhaseState p1 = *particle;
	MEGNO_Auxilary_Variables m1 = *megno;
	double deltaSq = p1.dx*p1.dx + p1.dy*p1.dy + p1.dvx*p1.dvx + p1.dvy*p1.dvy;
	double deltad_delta = p1.dx*p1.dvx + p1.dy*p1.dvy + p1.dax*p1.dvx + p1.day*p1.dvy;

	(*megno).Y = m1.Y +  deltad_delta / deltaSq * t * dt;
	(*megno).W = m1.W + dt * 2 * m1.Y / t ;
	(*megno).megno = m1.W/t;
}
void compute_variational_accs(PhaseState* restrict particle, double const mu1, double const mu2, double const Omega2, double t){
	const PhaseState p1 = *particle;
	// test-particle planet-1 separation
	const double dx1 = 1. - p1.x ;
	const double rho1sq = dx1 * dx1  + p1.y * p1.y ;
	const double rho1 = sqrt( rho1sq );
	const double rho1_inv3 = 1./(rho1sq * rho1);
	const double rho1_inv5 = rho1_inv3/rho1sq;
	// planet-2 location
	const double r2 = pow(Omega2, -2./3.);
	const double theta2 = Omega2 * t - t;
	const double x2 =  r2 * cos(theta2) ;
	const double y2 =  r2 * sin(theta2) ;	
	const double dx2 =  x2-p1.x ;
	const double dy2 =  y2-p1.y ;	
	const double rho2sq = dx2*dx2 + dy2*dy2;
	const double rho2 = sqrt(rho2sq);
	const double rho2_inv3 = 1./rho2sq / rho2;
	const double rho2_inv5 = rho2_inv3/rho2sq;
	// test-particle r	
	const double rsq = p1.x * p1.x + p1.y * p1.y;
	const double r = sqrt(rsq);
	const double r_inv3 = 1. / rsq / r; 
	const double r_inv5 = r_inv3 / rsq; 
	
	const double DFCENxx = 3 * p1.x * p1.x * r_inv5 - r_inv3;
	const double DFCENxy = 3 * p1.x * p1.y * r_inv5 ;
	const double DFCENyy = 3 * p1.y * p1.y * r_inv5 - r_inv3 ;
	const double DF1xx = 3.* dx1*dx1 * rho1_inv5 - rho1_inv3 \
					 	- 15.*p1.x*p1.x*p1.x * r_inv5 / rsq \
					 	+ 9*p1.x * r_inv5;
	const double DF1xy = -3. * dx1 * p1.y * rho1_inv5 \
						- 15.*p1.x * p1.x * p1.y * r_inv5 / rsq \
						+ 3*p1.y * r_inv5;
	const double DF1yy = ( 2*p1.y*p1.y - dx1*dx1 ) * rho1_inv5 + 3*p1.x*(p1.x*p1.x-4*p1.y*p1.y) * r_inv5 / rsq;
	const double DF2xx = 3*dx2*dx2*rho2_inv5-rho2_inv3;
	const double DF2xy = 3*dx2*dy2*rho2_inv5;
	const double DF2yy = 3*dy2*dy2*rho2_inv5-rho2_inv3;

	
	(*particle).dax =  p1.dvy + (DFCENxx + mu1*DF1xx+ mu2*DF2xx) * p1.dx + (DFCENxy + mu1*DF1xy+ mu2*DF2xy) * p1.dy ;
	(*particle).day = -p1.dvx + (DFCENxy + mu1*DF1xy+ mu2*DF2xy) * p1.dx + (DFCENyy + mu1*DF1yy+ mu2*DF2yy) * p1.dy ;
	
}

