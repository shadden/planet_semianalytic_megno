// Evolve under action of Hamiltonian
// H0 =  - 4/L^2 + n1 * (0.5 (X^2 + Y^2) - L)
// i.e. Kelerian Hamiltonian after transforming angle
// va
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "Resonances.h"
#include "ActionAngleIntegrate.h"


#define RT2 1.414213562373095
#define RT2INV 0.7071067811865475
#define INDX(ROW,COL) 4 * ROW + COL

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
	InitializeActionAngleSimulation(&sim,NresIn, innerResArr,NresOut,outerResArr,mu1,mu2,n1,n2,0.,0.,0.,L0,l0,X0,Y0,dt);
		
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
void InitializeActionAngleSimulation(ActionAngleSimulation* sim, int NresIn, int* innerRes,int NresOut, int* outerRes,\
									 double mu1,double mu2,double n1,double n2,double e1,double e2,double varpi2,\
									 double L0,double l0,double X0, double Y0,double dt)
									 {
	ActionAnglePhaseState state;
	MEGNO_Auxilary_Variables megno;
	ResonanceData r,r1;
	
	double alphaIn0=pow(1./n1,1./1.5);
	double alphaOut0=pow(n2,1./1.5);

	initialize_pars(&(sim->parameters),mu1,mu2,n1,n2,e1,e2,varpi2);
	ActionAnglePhaseStateInitialize(&(sim->state),  L0,  l0,  X0,  Y0);
	intialize_megno_vars(&(sim->megno));
	
	initialize_ResonanceData(&(sim->rIn));
	initialize_ResonanceData(&(sim->rOut));
	
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

//	printf("from C:\t dt: %.5f \t sim->dt: %.5f\n",dt,sim->dt);


}
void intialize_megno_vars( MEGNO_Auxilary_Variables* megno){

	(*megno).Y=0;
	(*megno).W=0;
	(*megno).megno=0;

}
void SimulationStep(ActionAngleSimulation* restrict sim){
		const double dt = sim->dt;
		const double t = sim->t;
		ActionAngle_H0_Advance( &(sim->state) , &(sim->parameters), 0.5 * (sim->dt));
		ActionAngle_H1_Advance_StormerVerlet(&(sim->state), &(sim->rIn) , &(sim->rOut) , &(sim->parameters), t ,dt);
		ActionAngle_H0_Advance( &(sim->state) , &(sim->parameters), 0.5 * dt);
		ActionAngle_Get_Var_Dot(&(sim->state), &(sim->rIn) , &(sim->rOut) , &(sim->parameters), t+dt);
		ActionAngle_update_megno_eqations(&(sim->state) , &(sim->megno) ,t+dt,dt);
		sim->t +=dt;
//		printf("from C:\t dt: %.5f \t sim->dt: %.5f \t sim->t: %.5f\n",dt,sim->dt,sim->t);

}
void initialize_pars(SimulationParameters* pars,double mu1,double mu2,double n1,double n2,double e1,double e2,double varpi2){
	pars->mu1=mu1;
	pars->mu2=mu2;
	pars->n1=n1;
	pars->n2=n2;
	pars->e1=e1;
	pars->e2=e2;
	pars->varpi2=varpi2;
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
void ActionAngle_H0_Advance( ActionAnglePhaseState* restrict Z ,SimulationParameters* restrict pars,const double dt){
		
	const ActionAnglePhaseState state = *Z;
	const double n1 = pars->n1;
	// Equations
	const double n1dt = n1 * dt;
	(*Z).l += dt * ( 8./(state.L*state.L*state.L) - n1 );
	(*Z).X  = cos(n1dt) * state.X - sin(n1dt) * state.Y;
	(*Z).Y  = sin(n1dt) * state.X + cos(n1dt) * state.Y;

	// Variationals
	(*Z).dl += dt *  -24./(state.L*state.L*state.L*state.L) * state.dL ;
	(*Z).dY  = sin(n1dt) * state.dX + cos(n1dt) * state.dY;
	(*Z).dX  = cos(n1dt) * state.dX - sin(n1dt) * state.dY;

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
	const double varpi2 = pars->varpi2;
	const double var[4] = {(Z->dl),(Z->dY),(Z->dL),(Z->dX)};

	double Xold	= Z->X;
	double Yold	= Z->Y;
	

	
	H1_Inner_Derivs(derivsIn,jacobianIn,Z,rIn ,mu1, n1,  e1, t);
	H1_Outer_Derivs(derivsOut,jacobianOut,Z,rOut,mu2,n1,n2,e2,varpi2,t);

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
	H1_Outer_Derivs(derivsOut,jacobianOut,Z,rOut,mu2,n1,n2,e2,varpi2,t);
	

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

//	Z->dldot += 0.5 * ( vardotOld[0] + vardotNew[0] );
//	Z->dYdot += 0.5 * ( vardotOld[1] + vardotNew[1] );
//	Z->dLdot += 0.5 * ( vardotOld[2] + vardotNew[2] );
//	Z->dXdot += 0.5 * ( vardotOld[3] + vardotNew[3] );

	
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
	const double varpi2 = pars->varpi2;

	ActionAnglePhaseState p1 = *state;


	(*state).dldot =  -24./(p1.L*p1.L*p1.L*p1.L) * p1.dL ;
	(*state).dYdot  = n1 * p1.dX ;
	(*state).dLdot  = 0. ;
	(*state).dXdot  = -n1 * p1.dY;
	H1_Inner_Derivs(derivsIn,jacobianIn,&p1,rIn ,mu1, n1,  e1, t);
	H1_Outer_Derivs(derivsOut,jacobianOut,&p1,rOut,mu2,n1,n2,e2,varpi2,t);
	
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
	
	double alpha = pow(n1,-2./3.);
	double rsq = 1 + alpha*alpha - 2 * alpha * cos(l0);
	double r = sqrt(rsq);
	Ldot += -2 * mu1 * (alpha * sin(l0) / rsq / r -  alpha * sin(l0)  );
	DLdotDl += -2 * mu1 * ( alpha * cos(l0) / rsq / r - 3 * alpha * alpha * sin(l0) * sin(l0) / rsq / rsq / r -  alpha * cos(l0) );
	
	// Add up inner planet effects 
	
	for(int i=0; i<NresIn; i++){

		j = *(rIn->ResonanceIndices + 3*i );
		o  = *(rIn->ResonanceIndices + 3*i + 1 );
		p = *(rIn->ResonanceIndices + 3*i + 2 );

		coeff = *( rIn->ResonanceCoefficients + ( MAX_ORDER + 1 )*i + p );

		theta = j * l0 + p * n1 * t + (o - p - 1) * g0;
		costheta = cos(theta);
		sintheta = sin(theta);

		// Derivatives
		factor  = o >= p+1 ? -RT2 * mu1 * coeff * (o-p) * mpow(e1,p) * mpow(E,o-p-1) : 0;		
		Ydot += factor * costheta;
		Xdot += factor * sintheta;
		Ldot +=  -2 * mu1  * coeff * j * mpow(e1,p) * mpow(E,o-p) * sin(theta + g0);

		// Variationals
		DYdotDl += -factor * j * sintheta;
		DXdotDl +=  factor * j * costheta;
		DLdotDl += -2 * mu1  * coeff * j*j * mpow(e1,p) * mpow(E,o-p) * cos(theta + g0);		

		factor  = o >= p+2 ? -2 * mu1 * coeff * (o-p) * (o-p-1) * mpow(e1,p) * mpow(E,o-p-2) : 0;

		DXdotDX +=  0.5 * factor * sin(theta-g0) ;	
		DYdotDX +=  0.5 * factor * cos(theta-g0) ;
		DXdotDY +=  0.5 * factor * cos(theta-g0) ;
		

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

void H1_Outer_Derivs(double* derivs,double* jacobian, ActionAnglePhaseState* Z, ResonanceData* restrict rOut ,const double mu2, const double n1,const double n2, const double e2,const double varpi2, const double t){

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

	double psi =  l0 - Dn2 * t  ;
	double rsq = 1 + alpha*alpha - 2 * alpha * cos(psi);
	double r = sqrt(rsq);

	Ldot += -2 * alpha * mu2 * ( alpha * sin(psi) / rsq / r -  alpha * sin(psi) );
	DLdotDl += -2 * alpha * mu2 * ( alpha * cos(psi) / rsq / r - 3 * alpha * alpha * sin(psi) * sin(psi) / rsq / rsq / r - alpha * cos(psi));



	// Add up resonances 
	for(int i=0; i<NresOut; i++){

		j = *(rOut->ResonanceIndices + 3*i );
		o  = *(rOut->ResonanceIndices + 3*i + 1 );
		p = *(rOut->ResonanceIndices + 3*i + 2 );

		coeff = *( rOut->ResonanceCoefficients + ( MAX_ORDER + 1 )*i + p );

		theta = j * Dn2 * t + (o-j) * l0 + (p-1) * g0 + (o-p) * (varpi2 + n1 * t);
		costheta=cos(theta);
		sintheta=sin(theta);
		
		factor = p>=1 ? -alpha * mu2 * RT2 * coeff * p * mpow(E,p-1) * mpow(e2,o-p) : 0 ;

		Ydot +=  factor * costheta;
		Xdot +=  factor * sintheta;
		Ldot +=  -2 * alpha * mu2  * coeff * (o-j) * mpow(e2,o-p) * mpow(E,p) * sin(theta + g0);
		
		// Variationals
		DYdotDl += -factor * (o-j) * sintheta;
		DXdotDl +=  factor * (o-j) * costheta;
		DLdotDl += -2 * alpha * mu2  * coeff * (o-j)* (o-j) * mpow(e2,o-p) * mpow(E,p) * cos(theta + g0);	
		
		factor = p>=2 ? -alpha * mu2 * 2 * coeff * p * (p-1) * mpow(E,p-2) * mpow(e2,o-p) : 0 ;

		DXdotDX +=  0.5 * factor * sin(theta-g0) ;	
		DYdotDX +=  0.5 * factor * cos(theta-g0) ;
		DXdotDY +=  0.5 * factor * cos(theta-g0) ;
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
