#include "shared.h"

// Fast inverse factorial lookup table
static const double invfactorial[35] = {1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000., 1./10888869450418352160768000000., 1./304888344611713860501504000000., 1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.};

static inline double fastabs(double x){
        return (x > 0.) ? x : -x;
}

static void stumpff_cs(double *restrict cs, double z) {
    unsigned int n = 0;
    while(fastabs(z)>0.1){
        z = z/4.;
        n++;
    }
    const int nmax = 15;
    double c_odd  = invfactorial[nmax];
    double c_even = invfactorial[nmax-1];
    for(int np=nmax-2;np>=5;np-=2){
        c_odd  = invfactorial[np]    - z *c_odd;
        c_even = invfactorial[np-1]  - z *c_even;
    }
    cs[5] = c_odd;
    cs[4] = c_even;
    cs[3] = invfactorial[3]  - z *cs[5];
    cs[2] = invfactorial[2]  - z *cs[4];
    cs[1] = invfactorial[1]  - z *cs[3];
    for (;n>0;n--){ 
        z = z*4.;
        cs[5] = (cs[5]+cs[4]+cs[3]*cs[2])*0.0625;
        cs[4] = (1.+cs[1])*cs[3]*0.125;
        cs[3] = 1./6.-z*cs[5];
        cs[2] = 0.5-z*cs[4];
        cs[1] = 1.-z*cs[3];
    }
    cs[0] = invfactorial[0]  - z *cs[2];
}

static void stumpff_cs3(double *restrict cs, double z) {
    unsigned int n = 0;
    while(fabs(z)>0.1){
        z = z/4.;
        n++;
    }
    const int nmax = 13;
    double c_odd  = invfactorial[nmax];
    double c_even = invfactorial[nmax-1];
    for(int np=nmax-2;np>=3;np-=2){
        c_odd  = invfactorial[np]    - z *c_odd;
        c_even = invfactorial[np-1]  - z *c_even;
    }
    cs[3] = c_odd;
    cs[2] = c_even;
    cs[1] = invfactorial[1]  - z *c_odd;
    cs[0] = invfactorial[0]  - z *c_even;
    for (;n>0;n--){ 
        cs[3] = (cs[2]+cs[0]*cs[3])*0.25;
        cs[2] = cs[1]*cs[1]*0.5;
        cs[1] = cs[0]*cs[1];
        cs[0] = 2.*cs[0]*cs[0]-1.;
    }
}

static void stiefel_Gs(double *restrict Gs, double beta, double X) {
    double X2 = X*X;
    stumpff_cs(Gs, beta*X2);
    Gs[1] *= X; 
    Gs[2] *= X2; 
    double _pow = X2*X;
    Gs[3] *= _pow; 
    _pow *= X;
    Gs[4] *= _pow; 
    _pow *= X;
    Gs[5] *= _pow; 
    return;
}

static void stiefel_Gs3(double *restrict Gs, double beta, double X) {
    double X2 = X*X;
    stumpff_cs3(Gs, beta*X2);
    Gs[1] *= X; 
    Gs[2] *= X2; 
    Gs[3] *= X2*X;
    return;
}

#define WHFAST_NMAX_QUART 64    ///< Maximum number of iterations for quartic solver
#define WHFAST_NMAX_NEWT  32    ///< Maximum number of iterations for Newton's method

void kepler_2D_advance(  PhaseState* restrict particle ,double _dt){
	
	const PhaseState p1 = *particle;
	const double r0 = sqrt(p1.x*p1.x + p1.y*p1.y);
	const double r0i = 1./r0;
	const double v2 =  p1.vx*p1.vx + p1.vy*p1.vy ;
	const double beta = 2.*r0i - v2;
	const double eta0 = p1.x*p1.vx + p1.y*p1.vy;
	const double zeta0 = 1 - beta*r0;
	double X;
	double Gs[6]; 
	double invperiod;  // only used for beta>0.
	double X_per_period = nan(""); // only used for beta>0. nan triggers Newton's method for beta<0.
	
	if (beta>0.){
        // Elliptic orbit
        const double sqrt_beta = sqrt(beta);
        invperiod = sqrt_beta*beta/(2.*M_PI);
        X_per_period = 2.*M_PI/sqrt_beta;
        const double dtr0i = _dt*r0i;
        X = dtr0i * (1. - dtr0i*eta0*0.5*r0i); // second order guess        
    }else{
        // Hyperbolic orbit
        X = 0.; // Initial guess 
    }
    
    unsigned int converged = 0;
    double oldX = X; 

    // Do one Newton step
    stiefel_Gs3(Gs, beta, X);
    const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
    double ri = 1./(r0 + eta0Gs1zeta0Gs2);
    X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+_dt);
	
	// Newton's method
	double oldX2 = nan("");             
	for (int n_hg=1;n_hg<WHFAST_NMAX_NEWT;n_hg++){
		oldX2 = oldX;
		oldX = X;
		stiefel_Gs3(Gs, beta, X);
		const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
		ri = 1./(r0 + eta0Gs1zeta0Gs2);
		X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+_dt);
	   
		if (X==oldX||X==oldX2){
			// Converged. Exit.
			converged = 1;
			break; 
		}
	}
    
        
    // If solver did not work, fallback to bisection 
    if (converged == 0){ 
        double X_min, X_max;
        if (beta>0.){
            //Elliptic
            X_min = X_per_period * floor(_dt*invperiod);
            X_max = X_min + X_per_period;
        }else{
            //Hyperbolic
            double h2 = r0*r0*v2-eta0*eta0;
            double q = h2/(1.+sqrt(1.-h2*beta));
            double vq = sqrt(h2)/q;
            X_min = 1./(vq+r0/_dt);
            X_max = _dt/q;
        }
        X = (X_max + X_min)/2.;
        do{
            stiefel_Gs3(Gs, beta, X);
            double s   = r0*X + eta0*Gs[2] + zeta0*Gs[3]-_dt;
            if (s>=0.){
                X_max = X;
            }else{
                X_min = X;
            }
            X = (X_max + X_min)/2.;
        }while (fastabs((X_max-X_min)/X_max)>1e-15);
        const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
        ri = 1./(r0 + eta0Gs1zeta0Gs2);
    }
    if (isnan(ri)){
        // Exception for (almost) straight line motion in hyperbolic case
        ri = 0.;
        Gs[1] = 0.;
        Gs[2] = 0.;
        Gs[3] = 0.;
    }


	// Note: These are not the traditional f and g functions.
    double f = -Gs[2]*r0i;
    double g = _dt - Gs[3];
    double fd = -Gs[1]*r0i*ri; 
    double gd = -Gs[2]*ri; 
        
    (*particle).x += f*p1.x + g*p1.vx;
    (*particle).y += f*p1.y + g*p1.vy;
        
    (*particle).vx += fd*p1.x + gd*p1.vx;
    (*particle).vy += fd*p1.y + gd*p1.vy;
    
    // Variations
    stiefel_Gs(Gs, beta, X);    // Recalculate (to get Gs[4] and Gs[5])
	
	double dr0 = (p1.x*p1.dx + p1.y*p1.dy)*r0i;
	double dbeta = -2.*dr0*r0i*r0i - 2.* (p1.vx*p1.dvx + p1.vy*p1.dvy);	
	double deta0 = p1.dx*p1.vx + p1.dy*p1.vy + p1.x*p1.dvx + p1.y*p1.dvy;
	double dzeta0 = -beta*dr0 - r0*dbeta;
	double G3beta = 0.5*(3.*Gs[5]-X*Gs[4]);
	double G2beta = 0.5*(2.*Gs[4]-X*Gs[3]);
	double G1beta = 0.5*(Gs[3]-X*Gs[2]);
	double tbeta = eta0*G2beta + zeta0*G3beta;
	double dX = -1.*ri*(X*dr0 + Gs[2]*deta0+Gs[3]*dzeta0+tbeta*dbeta);
	double dG1 = Gs[0]*dX + G1beta*dbeta; 
	double dG2 = Gs[1]*dX + G2beta*dbeta;
	double dG3 = Gs[2]*dX + G3beta*dbeta;
	double dr = dr0 + Gs[1]*deta0 + Gs[2]*dzeta0 + eta0*dG1 + zeta0*dG2;
	double df = Gs[2]*dr0*r0i*r0i - dG2*r0i;
	double dg = -dG3;
	double dfd = -dG1*r0i*ri + Gs[1]*(dr0*r0i+dr*ri)*r0i*ri;
	double dgd = -dG2*ri + Gs[2]*dr*ri*ri;

	(*particle).dx += f*p1.dx + g*p1.dvx + df*p1.x + dg*p1.vx;
	(*particle).dy += f*p1.dy + g*p1.dvy + df*p1.y + dg*p1.vy;

	(*particle).dvx += fd*p1.dx + gd*p1.dvx + dfd*p1.x + dgd*p1.vx;
	(*particle).dvy += fd*p1.dy + gd*p1.dvy + dfd*p1.y + dgd*p1.vy;
}
