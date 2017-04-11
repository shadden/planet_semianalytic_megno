// Evolve under action of Hamiltonian
// H0 =  - 4/L^2 + n1 * (0.5 (X^2 + Y^2) - L)
// i.e. Kelerian Hamiltonian after transforming angle
// va
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "shared.h"
#include "LaplaceCoefficients.h"
#include "ActionAngleIntegrate.h"


#define RT2 1.414213562373095
#define RT2INV 0.7071067811865475
#define INDX(ROW,COL) 4 * ROW + COL
#define PRINT 0


double const symplecticJ[4][4] = {{  0,  0,  1,  0},\
								  {  0,  0,  0,  1},\
								  { -1,  0,  0,  0},\
								  {  0, -1,  0,  0}};
double mpow (double, int);
double mpow2(double, int);
void freeSimluation(ActionAngleSimulation* sim);

int IntegrateSimulation(ActionAngleSimulation* sim, const double tFin){
	const double t0 = sim->t;
	const double dt = sim->dt;
	const int Nstep = (int)((tFin - t0)/dt);

	for(int i=0; i<Nstep; i++){	
		SimulationStep(sim);
	}
	sim->t = t0 + Nstep * dt;
	return Nstep;
}

double CircularFirstOrderResonanceMEGNOIntegration
(int NresIn, int* innerRes, int NresOut, int* outerRes,double mu1,double mu2,double n1,double n2, double tFin)
 {
	double L0=2.,l0=0,X0=0,Y0=0;
	ActionAngleSimulation sim;
	const double dt = 2*M_PI / 20. ;
	int innerResArr[3*NresIn];
	int outerResArr[3*NresOut];

	for(int j=0;j<NresIn;j++){
		innerResArr[3*j]=*(innerRes+j);
		innerResArr[3*j+1]=1;
		innerResArr[3*j+2]=0;
	}
	for(int j=0;j<NresOut;j++){
		outerResArr[3*j]=*(outerRes+j);
		outerResArr[3*j+1]=1;
		outerResArr[3*j+2]=1;
	}
	InitializeActionAngleSimulation(&sim,NresIn,true, innerResArr,NresOut,true,outerResArr,mu1,mu2,n1,n2,0.,0.,0.,0.,L0,l0,X0,Y0,dt);
		
	int Nstep;
	Nstep = (int)(tFin/dt);
	for(int i=0; i<Nstep; i++){	
		SimulationStep(&sim);
	}
	double megno = sim.megno.megno;	
	freeSimluation(&sim);
	return megno;
 }
 
 void freeSimluation(ActionAngleSimulation* sim){
 	free_ResonanceData(&(sim->rIn));
 	 free_ResonanceData(&(sim->rOut));
 }
 
static inline double random_real(){
        return -1+2.*((float)rand())/RAND_MAX;
}
void InitializeActionAngleSimulation(ActionAngleSimulation* sim, int NresIn, int includeInnerZeroth, int* innerRes,int NresOut,int includeOuterZeroth, int* outerRes,\
									 double mu1,double mu2,double n1,double n2,double e1,double e2,double lambda2, double varpi2,\
									 double L0,double l0,double X0, double Y0,double dt)
									 {
	ActionAnglePhaseState state;
	MEGNO_Auxilary_Variables megno;
	ResonanceData r,r1;
	
	double alphaIn0=pow(1./n1,1./1.5);
	double alphaOut0=pow(n2,1./1.5);

	initialize_pars(&(sim->parameters),mu1,mu2,n1,n2,e1,e2,lambda2, varpi2);
	ActionAnglePhaseStateInitialize(&(sim->state),  L0,  l0,  X0,  Y0);
	intialize_megno_vars(&(sim->megno));
	
	initialize_ResonanceData(&(sim->rIn),includeInnerZeroth);
	initialize_ResonanceData(&(sim->rOut),includeOuterZeroth);
	
	int j,o,p;
	for(int i=0;i<NresIn;i++){
		j = *(innerRes + 3*i + 0 );
		o = *(innerRes + 3*i + 1 );
		p = *(innerRes + 3*i + 2 );
		AddResonance(&(sim->rIn) ,j,o,p, alphaIn0 );
	}
	for(int i=0;i<NresOut;i++){
		j = *(outerRes + 3*i + 0 );
		o = *(outerRes + 3*i + 1 );
		p = *(outerRes + 3*i + 2 );
		AddResonance(&(sim->rOut) ,j,o,p, alphaOut0 );
	}
	
	sim->t=0;
	sim->dt = dt;



}
void intialize_megno_vars( MEGNO_Auxilary_Variables* megno){

	(*megno).Y=0;
	(*megno).W=0;
	(*megno).megno=0;

}
void SimulationStep(ActionAngleSimulation* restrict sim){
		const double dt = sim->dt;
		const double t = sim->t;
		ActionAngle_H0_Advance( &(sim->state) , &(sim->parameters), t, 0.5 * (sim->dt));
		ActionAngle_H1_Advance_StormerVerlet(&(sim->state), &(sim->rIn) , &(sim->rOut) , &(sim->parameters), t ,dt);
		ActionAngle_H0_Advance( &(sim->state) , &(sim->parameters), t + 0.5 * dt ,0.5 * dt);
		ActionAngle_Get_Var_Dot(&(sim->state), &(sim->rIn) , &(sim->rOut) , &(sim->parameters), t+dt);
		ActionAngle_update_megno_eqations(&(sim->state) , &(sim->megno) ,t+dt,dt);
		sim->t +=dt;

}


void initialize_pars(SimulationParameters* pars,double mu1,double mu2,double n1,double n2,double e1,double e2,double lambda2,double varpi2){
	pars->mu1=mu1;
	pars->mu2=mu2;
	pars->n1=n1;
	pars->n2=n2;
	pars->e1=e1;
	pars->e2=e2;
	pars->lambda2=lambda2;
	pars->varpi2=varpi2;
	
	pars->alphaIn=pow(n1,-2./3.);
 	pars->alphaOut = pow(n2,2./3.);
 	
 	// secular coefficients
 	pars->fSecIn = secularF2(pars->alphaIn);
 	pars->gSecIn = secularF10(pars->alphaIn);
 	pars->fSecOut = secularF2(pars->alphaOut);
 	pars->gSecOut = secularF10(pars->alphaOut);

}
void ActionAnglePhaseStateInitialize(ActionAnglePhaseState* restrict Z, double L0, double l0, double X0, double Y0){
	Z->L = L0;
	Z->l = l0;
	Z->X = X0;
	Z->Y = Y0;
	double dvec[4];
	double normsq=0;
	
	// generate random reals	
	for(int i=0; i < 4; i++){
		dvec[i] = random_real();
		normsq += dvec[i]*dvec[i];
	}
	const double EPS = 1.e-14;
	(*Z).dL = EPS*dvec[0]/sqrt(normsq);
	(*Z).dl = EPS*dvec[1]/sqrt(normsq);
	(*Z).dX = EPS*dvec[2]/sqrt(normsq);
	(*Z).dY = EPS*dvec[3]/sqrt(normsq);
	(*Z).dLdot = 0 ;
	(*Z).dldot = 0 ;
	(*Z).dXdot = 0 ;
	(*Z).dYdot = 0 ;


}

void ActionAngle_H0_Advance( ActionAnglePhaseState* restrict Z ,SimulationParameters* restrict pars, const double t,const double dt){
		
	const ActionAnglePhaseState state = *Z;
	const double n1 = pars->n1;
	double const varpi2 = pars->varpi2;

	double const fIn =pars->fSecIn;
	double const gIn =pars->gSecIn;

	double const fOut =pars->fSecOut;
	double const gOut =pars->gSecOut;
	
	double const mu1 =pars->mu1;
	double const mu2 =pars->mu2;
	double const alpha2 =pars->alphaOut;
	
	double const e1 =pars->e1;
	double const e2 =pars->e2;

	// Equations

	
	double X_free, Y_free, x_forced, y_forced, w_sec,X_forced, Y_forced;
	double x1 = sqrt(2.) * e1 * cos(n1*t);
	double y1 = sqrt(2.) * e1 * sin(n1*t);
	double x2 = sqrt(2.) * e2 * cos(n1 * t - varpi2);
	double y2 = sqrt(2.) * e2 * sin(n1 * t - varpi2);
	
	
	// NOTE: this is gamma_dot = -1 * pomega_dot
	w_sec = -2 * alpha2 * mu2 * fOut - 2 * mu1 * fIn;

	const double dtheta = (n1+w_sec) * dt;
	const double n1dt =  (n1) * (dt);
	const double Sn1dt =  sin(n1dt);
	const double Cn1dt =  cos(n1dt);
	const double Sdtheta = sin(dtheta);
	const double Cdtheta = cos(dtheta);
	
	x_forced =  ( mu1 * gIn * x1  + alpha2 * mu2 * gOut * x2) / w_sec;
	y_forced =  ( mu1 * gIn * y1  + alpha2 * mu2 * gOut * y2) / w_sec;


	
	X_free = state.X - x_forced;
	Y_free = state.Y - y_forced;
	


	(*Z).l += dt * ( 8./(state.L*state.L*state.L) - n1 );
	(*Z).X  = Cdtheta * X_free - Sdtheta * Y_free + x_forced * Cn1dt - y_forced * Sn1dt;
	(*Z).Y  = Sdtheta * X_free + Cdtheta * Y_free + x_forced * Sn1dt + y_forced * Cn1dt;
	

	// Variationals
	(*Z).dl += dt *  -24./(state.L*state.L*state.L*state.L) * state.dL ;
	(*Z).dY  = Sdtheta * state.dX + Cdtheta * state.dY;
	(*Z).dX  = Cdtheta * state.dX - Sdtheta * state.dY;

}
void ActionAngle_H1_Advance_StormerVerlet(ActionAnglePhaseState* Z, ResonanceData* restrict rIn, ResonanceData* restrict rOut ,SimulationParameters* restrict pars, const double t, const double dt){
	
	double* derivsIn = (double *)malloc(4 * sizeof(double));
	double* jacobianIn = (double *)malloc(4*4 * sizeof(double));
	double* derivsOut = (double *)malloc(4 * sizeof(double));
	double* jacobianOut = (double *)malloc(4*4 * sizeof(double));
	double dt_2 = 0.5 * dt;

	double Xdot,XdotOld;
	double Ydot,YdotOld;
	double Ldot ;
	
	
	const double mu1 = pars->mu1;
	const double mu2 = pars->mu2;
	const double n1 = pars->n1;
	const double n2 = pars->n2;
	const double e1 = pars->e1;
	const double e2 = pars->e2;
	const double lambda2 = pars->lambda2;
	const double varpi2 = pars->varpi2;
	const double var[4] = {(Z->dl),(Z->dY),(Z->dL),(Z->dX)};

	double Xold	= Z->X;
	double Yold	= Z->Y;
	

	
	H1_Inner_Derivs(derivsIn,jacobianIn,Z,rIn ,mu1, n1,  e1, t);
	H1_Outer_Derivs(derivsOut,jacobianOut,Z,rOut,mu2,n1,n2,e2,lambda2,varpi2,t);

	YdotOld = *(derivsIn+1) + *(derivsOut+1); 
	XdotOld = *(derivsIn+3) + *(derivsOut+3);

	Z->Y += dt * YdotOld;
	Z->X += dt * XdotOld;
	Ldot = *(derivsIn+2) + *(derivsOut+2);
	Z->L += dt_2 * Ldot;
	

	double vardotOld[4];
	for(int j=0; j<4;j++)
	{
		vardotOld[j]=0;
		for(int k=0; k<4;k++){
			 vardotOld[j] += (*(jacobianIn + INDX(j,k)) +  *(jacobianOut + INDX(j,k)) )  * var[k];
		}
	}
	
	H1_Inner_Derivs(derivsIn,jacobianIn,Z,rIn ,mu1, n1,  e1, t);
	H1_Outer_Derivs(derivsOut,jacobianOut,Z,rOut,mu2,n1,n2,e2,lambda2,varpi2,t);
	

	Ldot = *(derivsIn+2) + *(derivsOut+2) ;
	Ydot = *(derivsIn+1) + *(derivsOut+1) ; 
	Xdot = *(derivsIn+3) + *(derivsOut+3) ;
	double vardotNew[4];
	for(int j=0; j<4;j++)
	{
		vardotNew[j]=0;
		for(int k=0; k<4;k++){
			 vardotNew[j] += (*(jacobianIn + INDX(j,k)) +  *(jacobianOut + INDX(j,k)) )  * (var[k] + vardotOld[k] * dt);
		}
	}

	
	Z->L += dt_2 * Ldot;
	Z->X = Xold + dt_2 * (XdotOld + Xdot);
	Z->Y = Yold + dt_2 * (YdotOld + Ydot);
	
	Z->dl+= dt_2 * ( vardotOld[0] + vardotNew[0] );
	Z->dY+= dt_2 * ( vardotOld[1] + vardotNew[1] );
	Z->dL+= dt_2 * ( vardotOld[2] + vardotNew[2] );
	Z->dX+= dt_2 * ( vardotOld[3] + vardotNew[3] );
	
	free(derivsIn);
	free(derivsOut);
	free(jacobianIn);
	free(jacobianOut);

}
void ActionAngle_Get_Var_Dot(ActionAnglePhaseState* state, ResonanceData* restrict rIn, ResonanceData* restrict rOut ,SimulationParameters* restrict pars, const double t){
	
	double* derivsIn = (double *)malloc(4 * sizeof(double));
	double* jacobianIn = (double *)malloc(4*4 * sizeof(double));
	double* derivsOut = (double *)malloc(4 * sizeof(double));
	double* jacobianOut = (double *)malloc(4*4 * sizeof(double));	
	const double mu1 = pars->mu1;
	const double mu2 = pars->mu2;
	const double n1 = pars->n1;
	const double n2 = pars->n2;
	const double e1 = pars->e1;
	const double e2 = pars->e2;
	const double lambda2 = pars->lambda2;
	const double varpi2 = pars->varpi2;

	ActionAnglePhaseState p1 = *state;


	(*state).dldot =  -24./(p1.L*p1.L*p1.L*p1.L) * p1.dL ;
	(*state).dYdot  = n1 * p1.dX ;
	(*state).dLdot  = 0. ;
	(*state).dXdot  = -n1 * p1.dY;
	H1_Inner_Derivs(derivsIn,jacobianIn,&p1,rIn ,mu1, n1,  e1, t);
	H1_Outer_Derivs(derivsOut,jacobianOut,&p1,rOut,mu2,n1,n2,e2,lambda2,varpi2,t);
	
	double vardot[4];
	const double var[4] = {(p1.dl),(p1.dY),(p1.dL),(p1.dX)};
	
	for(int j=0; j<4;j++)
	{
		vardot[j]=0;
		for(int k=0; k<4;k++){
			 vardot[j] += (*(jacobianIn + INDX(j,k)) +  *(jacobianOut + INDX(j,k)) )  * var[k];
		}
	}
	(*state).dldot += vardot[0] ;
	(*state).dYdot += vardot[1] ;
	(*state).dLdot += vardot[2] ;
	(*state).dXdot += vardot[3] ;

	free(derivsIn);
	free(derivsOut);
	free(jacobianIn);
	free(jacobianOut);



}
void ActionAngle_update_megno_eqations(ActionAnglePhaseState* restrict particle , MEGNO_Auxilary_Variables* megno ,double t,double dt){
	const ActionAnglePhaseState p1 = *particle;
	MEGNO_Auxilary_Variables m1 = *megno;
	double deltaSq = p1.dl*p1.dl + p1.dY*p1.dY + p1.dL*p1.dL + p1.dX*p1.dX;
	double deltad_delta = p1.dl*p1.dldot + p1.dY*p1.dYdot + p1.dL*p1.dLdot + p1.dX*p1.dXdot;
	
	(*megno).Y = m1.Y +  deltad_delta / deltaSq * t * dt;
	(*megno).W = m1.W + dt * 2 * m1.Y / t ;
	(*megno).megno = m1.W/t;
}

void cosine_sine_array(const double cosx,const double sinx, const int Nmax, double* cosarray, double* sinarray){
	assert(Nmax>=0);
	*(cosarray) = 1;
	*(sinarray) = 0;
	for(int k=1; k<Nmax; k++){
		*(cosarray+k)  = cosx * (*(cosarray+k-1)) - sinx * (*(sinarray+k-1));
		*(sinarray+k)  = sinx * (*(cosarray+k-1)) + cosx * (*(sinarray+k-1));
	}
}

void H1_Inner_Derivs(double* derivs,double* jacobian ,ActionAnglePhaseState* Z, ResonanceData* restrict rIn ,const double mu1, const double n1, const double e1, const double t){

	const int NresIn = rIn->Nres;

	int j,o,p;
	double coeff;
	const double X0 = Z->X;
	const double Y0 = Z->Y;
	const double l0 = Z->l;

	const double Esq = 0.5*(X0*X0 + Y0*Y0);
	const double E = sqrt(Esq);
	const double g0 = atan2(Y0,X0);

	double Xdot=0;
	double Ydot=0;
	double Ldot=0;
	double theta,factor;
	double costheta,sintheta;
	double DLdotDl=0,DYdotDl=0,DXdotDl=0,DXdotDY=0,DYdotDX=0,DXdotDX=0;

	// Zero-th order resonance effects
	if (rIn->IncludeZeroth){
		double alpha = pow(n1,-2./3.);
		double c_l0 = cos(l0);
		double s_l0 = sin(l0);
		double rsq = 1 + alpha*alpha - 2 * alpha * c_l0;
		double r = sqrt(rsq);
		Ldot += -2 * mu1 * (alpha * s_l0 / rsq / r -  alpha * s_l0  );
		DLdotDl += -2 * mu1 * ( alpha * c_l0 / rsq / r - 3 * alpha * alpha * s_l0 * s_l0 / rsq / rsq / r -  alpha * c_l0 );
	}

	// Add up inner planet effects 	
	// get cosine and sine data for

	double cos_l_array[MAX_J]; 
	double sin_l_array[MAX_J]; 
	cosine_sine_array(cos(l0),sin(l0),(rIn->MaxJ)+1,cos_l_array,sin_l_array);

	double cos_g_array[MAX_ORDER+1]; 
	double sin_g_array[MAX_ORDER+1]; 
	cosine_sine_array(cos(g0),sin(g0),(rIn->MaxOrder)+1, cos_g_array,sin_g_array);
	
	double cos_n1t_array[MAX_ORDER+1]; 
	double sin_n1t_array[MAX_ORDER+1]; 
	cosine_sine_array(cos(n1*t),sin(n1*t),(rIn->MaxOrder)+1, cos_n1t_array,sin_n1t_array);

	double c_j_l;
	double s_j_l;
	double c_p_n1t;
	double s_p_n1t;
	double c_op1_g;
	double s_op1_g;
	double c_g = cos_g_array[1];
	double s_g = sin_g_array[1];
	
	for(int i=0; i<NresIn; i++){

		j = *(rIn->ResonanceIndices + 3*i );
		o  = *(rIn->ResonanceIndices + 3*i + 1 );
		p = *(rIn->ResonanceIndices + 3*i + 2 );

		coeff = *( rIn->ResonanceCoefficients + ( MAX_ORDER + 1 )*i + p );

		
		c_j_l = cos_l_array[j];
		s_j_l = sin_l_array[j];
		c_p_n1t = cos_n1t_array[p];
		s_p_n1t = sin_n1t_array[p];
		c_op1_g = o-p-1 >= 0 ? cos_g_array[o-p-1] : c_g;
		s_op1_g = o-p-1 >= 0 ? sin_g_array[o-p-1] : -1*s_g;
		theta = j * l0 + p * n1 * t + (o - p - 1) * g0;			
		costheta = c_j_l * ( c_p_n1t * c_op1_g - s_p_n1t * s_op1_g ) - s_j_l * ( c_p_n1t * s_op1_g + s_p_n1t * c_op1_g);		
		sintheta = c_j_l * ( c_p_n1t * s_op1_g + s_p_n1t * c_op1_g ) + s_j_l * ( c_p_n1t * c_op1_g - s_p_n1t * s_op1_g);
	
#if PRINT
		printf("inner %d %d %d: %g \t %g \n",j,o,p,costheta-cos(theta),sintheta-sin(theta) );
#endif
	
		// Derivatives
		factor  = o >= p+1 ? -RT2 * mu1 * coeff * (o-p) * mpow(e1,p) * mpow(E,o-p-1) : 0;		
		Ydot += factor * costheta;
		Xdot += factor * sintheta;
		Ldot +=  -2 * mu1  * coeff * j * mpow(e1,p) * mpow(E,o-p) * (sintheta * c_g + costheta * s_g );
		
		// Variationals
		DYdotDl += -factor * j * sintheta;
		DXdotDl +=  factor * j * costheta;
		DLdotDl += -2 * mu1  * coeff * j*j * mpow(e1,p) * mpow(E,o-p) * ( costheta*c_g-sintheta*s_g );		

		factor  = o >= p+2 ? -2 * mu1 * coeff * (o-p) * (o-p-1) * mpow(e1,p) * mpow(E,o-p-2) : 0;
		
		DXdotDX +=  0.5 * factor * (sintheta*c_g-costheta*s_g) ;	
 		DYdotDX +=  0.5 * factor * ((costheta * c_g) + (sintheta * s_g)) ;
 		DXdotDY +=  0.5 * factor * ((costheta * c_g) + (sintheta * s_g)) ;
		

	}
	// Coordinates z_i are:
	//	i	coord
	//	-	-----
	//	0	l
	//	1	Y
	//	2	L
	//	3	X
	
	*(derivs) = 0.;
	*(derivs+1) = Ydot;
	*(derivs+2) = Ldot;
	*(derivs+3) = Xdot;
	

	double Hij[4][4];
	Hij[0][0] =  -DLdotDl ; // H_l,l
	Hij[0][1] = -DXdotDl ; // H_l,Y
	Hij[0][2] =  0.	   ; // H_l,L		
	Hij[0][3] =  DYdotDl ; // H_l,X	

	Hij[1][1] = -DXdotDY ; // H_Y,Y
	Hij[1][2] = 0.       ; // H_Y,L
	Hij[1][3] = -DXdotDX  ; // H_Y,X

	Hij[2][2] = 0.       ; // H_L,L
	Hij[2][3] = 0.       ; // H_L,X

	Hij[3][3] = DYdotDX  ; // H_X,X
	
	int row,col;
	
	double jacobian_ij;
	for(row=0;row<4;row++){
		for(col=0;col<4;col++){
		jacobian_ij = 0;
		for(int l=0; l<4; l++){
			if(l>col){
				Hij[l][col] = Hij[col][l];
			}
			jacobian_ij+= symplecticJ[row][l] * Hij[l][col];
		}
				*( jacobian + INDX(row,col) ) = jacobian_ij;
		}
	}

}

void H1_Outer_Derivs(double* derivs,double* jacobian, ActionAnglePhaseState* Z, ResonanceData* restrict rOut ,const double mu2, const double n1,const double n2, const double e2,const double lambda2,const double varpi2, const double t){

	const int NresOut = rOut->Nres;

	int j,o,p;
	double coeff;
	
	const double Dn2 = n2 - n1;
	const double alpha =pow(n2, 1./1.5);
	
	const double X0 = Z->X;
	const double Y0 = Z->Y;
	const double l0 = Z->l;

	const double Esq = 0.5*(X0*X0 + Y0*Y0);
	const double E = sqrt(Esq);
	const double g0 = atan2(Y0,X0);
		

	double Xdot=0;
	double Ydot=0;
	double Ldot=0;
	double theta,factor;
	double costheta,sintheta;
	double DLdotDl=0,DYdotDl=0,DXdotDl=0,DXdotDY=0,DYdotDX=0,DXdotDX=0;
	
	// Zero-th order resonance effects

	if (rOut->IncludeZeroth){
		double psi =  l0 - Dn2 * t - lambda2 ;
		
		double c_psi = cos(psi);
		double s_psi = sin(psi);
		
		double rsq = 1 + alpha*alpha - 2 * alpha * c_psi;
		double r = sqrt(rsq);
	
		Ldot += -2 * alpha * mu2 * ( alpha * s_psi / rsq / r -  alpha * s_psi );
		DLdotDl += -2 * alpha * mu2 * ( alpha * c_psi / rsq / r - 3 * alpha * alpha * s_psi * s_psi / rsq / rsq / r - alpha * c_psi);
	}
	
	// get cosine and sine data

	double cos_l2_array[MAX_J]; 
	double sin_l2_array[MAX_J]; 
	cosine_sine_array(cos(Dn2 * t + lambda2),sin(Dn2 * t + lambda2),(rOut->MaxJ)+1,cos_l2_array,sin_l2_array);
	double cos_l_array[MAX_J]; 
	double sin_l_array[MAX_J]; 
	cosine_sine_array(cos(l0),sin(l0),(rOut->MaxJ)+1,cos_l_array,sin_l_array);


	double cos_g_array[MAX_ORDER+1]; 
	double sin_g_array[MAX_ORDER+1]; 
	cosine_sine_array(cos(g0),sin(g0),(rOut->MaxOrder)+1, cos_g_array,sin_g_array);
	
	double cos_g2_array[MAX_ORDER+1]; 
	double sin_g2_array[MAX_ORDER+1]; 
	cosine_sine_array(cos(n1 * t - varpi2),sin(n1 * t - varpi2),(rOut->MaxOrder)+1, cos_g2_array,sin_g2_array);

	double c_j_l2;
	double s_j_l2;
	double c_j_l0;
	double s_j_l0;
	double c_j_g;
	double s_j_g;
	double c_j_g2;
	double s_j_g2;
 	double c_g = cos_g_array[1];
 	double s_g = sin_g_array[1];
 	
	// Add up resonances
	for(int i=0; i<NresOut; i++){

		j = *(rOut->ResonanceIndices + 3*i );
		o  = *(rOut->ResonanceIndices + 3*i + 1 );
		p = *(rOut->ResonanceIndices + 3*i + 2 );

		coeff = *( rOut->ResonanceCoefficients + ( MAX_ORDER + 1 )*i + p );

		theta = j * (Dn2 * t + lambda2) + (o-j) * l0 + (p-1) * g0 + (o-p) * (n1 * t - varpi2);
		c_j_l2 = cos_l2_array[j];
		s_j_l2 = sin_l2_array[j];
		c_j_l0 = cos_l_array[j-o];
		s_j_l0 = -1*sin_l_array[j-o];
		c_j_g = p-1 >= 0 ? cos_g_array[p-1] : cos_g_array[1];
		s_j_g = p-1 >= 0 ? sin_g_array[p-1] : -1*sin_g_array[1];
		c_j_g2 =  cos_g2_array[o-p];
		s_j_g2 =  sin_g2_array[o-p];

		
		costheta = c_j_g*c_j_g2*c_j_l0*c_j_l2 - c_j_l0*c_j_l2*s_j_g*s_j_g2 - \
		c_j_g2*c_j_l2*s_j_g*s_j_l0 - c_j_g*c_j_l2*s_j_g2*s_j_l0 - \
		c_j_g2*c_j_l0*s_j_g*s_j_l2 - c_j_g*c_j_l0*s_j_g2*s_j_l2 - \
		c_j_g*c_j_g2*s_j_l0*s_j_l2 + s_j_g*s_j_g2*s_j_l0*s_j_l2;
		
		sintheta = c_j_g2*c_j_l0*c_j_l2*s_j_g + c_j_g*c_j_l0*c_j_l2*s_j_g2 + \
		c_j_g*c_j_g2*c_j_l2*s_j_l0 - c_j_l2*s_j_g*s_j_g2*s_j_l0 + \
		c_j_g*c_j_g2*c_j_l0*s_j_l2 - c_j_l0*s_j_g*s_j_g2*s_j_l2 - \
		c_j_g2*s_j_g*s_j_l0*s_j_l2 - c_j_g*s_j_g2*s_j_l0*s_j_l2 ;
#if PRINT
		printf("outer %d %d %d: %g \t %g \n",j,o,p,costheta-cos(theta),sintheta-sin(theta) );
#endif		
		factor = p>=1 ? -alpha * mu2 * RT2 * coeff * p * mpow(E,p-1) * mpow(e2,o-p) : 0 ;

		Ydot +=  factor * costheta;
		Xdot +=  factor * sintheta;
		Ldot +=  -2 * alpha * mu2  * coeff * (o-j) * mpow(e2,o-p) * mpow(E,p) * (sintheta*c_g + costheta*s_g);
		
		// Variationals
		DYdotDl += -factor * (o-j) * sintheta;
		DXdotDl +=  factor * (o-j) * costheta;
		DLdotDl += -2 * alpha * mu2  * coeff * (o-j)* (o-j) * mpow(e2,o-p) * mpow(E,p) * (costheta*c_g - sintheta*s_g);	

		factor = p>=2 ? -alpha * mu2 * 2 * coeff * p * (p-1) * mpow(E,p-2) * mpow(e2,o-p) : 0 ;

		DXdotDX +=  0.5 * factor * (sintheta*c_g - costheta*s_g);	
		DYdotDX +=  0.5 * factor * (costheta*c_g + sintheta*s_g) ;
		DXdotDY +=  0.5 * factor * (costheta*c_g + sintheta*s_g) ;
	}
	// Coordinates z_i are:
	//	i	coord
	//	-	-----
	//	0	l
	//	1	Y
	//	2	L
	//	3	X
	
	*(derivs) = 0.;
	*(derivs+1) = Ydot;
	*(derivs+2) = Ldot;
	*(derivs+3) = Xdot;
	

	double Hij[4][4];
	Hij[0][0] =  -DLdotDl ; // H_l,l
	Hij[0][1] = -DXdotDl ; // H_l,Y
	Hij[0][2] =  0.	   ; // H_l,L		
	Hij[0][3] =  DYdotDl ; // H_l,X	

	Hij[1][1] = -DXdotDY ; // H_Y,Y
	Hij[1][2] = 0.       ; // H_Y,L
	Hij[1][3] = -DXdotDX  ; // H_Y,X

	Hij[2][2] = 0.       ; // H_L,L
	Hij[2][3] = 0.       ; // H_L,X

	Hij[3][3] = DYdotDX  ; // H_X,X


	double jacobian_ij;
	for(int row=0;row<4;row++){
		for(int col=0;col<4;col++){
		jacobian_ij = 0;
		for(int l=0; l<4; l++){
			if(l>col){
				Hij[l][col] = Hij[col][l];
			}
			jacobian_ij+= symplecticJ[row][l] * Hij[l][col];
		}
				*( jacobian + INDX(row,col) ) = jacobian_ij;
		}
	}
	
}


double mpow(double b, int exp){
   if (exp < 0)
   {
      b = 1/b;
      exp *= -1;
   }
   
   return mpow2(b,exp);
}
double mpow2(double b, int exp){
   if (exp > 2)
   {
      if (exp%2)
         return b*mpow(b*b,(int)exp/2);
      else
         return mpow(b*b,exp/2);
   }
   else if (2 == exp)
      return b*b;
   else if (1 == exp)
      return b;

   return 1.0; // exp == 0
}
