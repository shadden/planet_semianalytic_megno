#include "shared.h"


int main(void){

MEGNO_Auxilary_Variables megno;
intialize_megno_vars( &megno);

PhaseState planet;
intialize_particle(&planet,1.05,0.0);

double mu1=1.0e-5;
double mu2=0.0;
double Omega2 = 0.51;
const double dt=2.*M_PI / 30. ;

for(int i=0; i< 1000*30; i++){

// Both planets 
	two_circular_perturbers_advance( mu1 , mu2 , Omega2 , &planet , i*dt , dt*0.5 );	
	kepler_2D_advance(&planet, dt);
	rotaion_advance(&planet,dt);
	two_circular_perturbers_advance( mu1 , mu2 , Omega2 , &planet , i*dt + dt , dt*0.5 );

	compute_variational_accs(&planet,mu1,mu2,Omega2,i*dt + dt);
	update_megno_eqations( &planet , &megno ,i*dt + dt, dt);

	if(i%50==0){
	printf("%.8f\t%.8g\t%.8g\t%.8g\n",(i)*dt,megno.Y,megno.W,megno.megno);
	}


	


}


return 0;
}
