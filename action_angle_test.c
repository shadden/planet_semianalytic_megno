#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "Resonances.h"
#include "ActionAngleIntegrate.h"
//#include "shared.h"


#define INTEGRATE 0
#define DERIVS 1
#define INDX(ROW,COL) 4 * ROW + COL

const double L0 = 2.;
const double l0= -1./5.;
const double X0=1;
const double Y0=1;

const double mu1 = 1;
const double mu2 = 1;
const double n1 = 7./5.;
const double n2 = 5./7.;
const double e1=0.;
const double e2=0.;
const double varpi2=0.;

double dt = 2*M_PI / 20.;
double tFin = 2*M_PI*1E4;

int Nstep; 


int main(){


	SimulationParameters pars;
	ActionAnglePhaseState state;
	MEGNO_Auxilary_Variables megno;
	ResonanceData r,r1;

	initialize_pars(&pars,mu1,mu2,n1,n2,e1,e2,varpi2);
	ActionAnglePhaseStateInitialize(&state,  L0,  l0,  X0,  Y0);
	intialize_megno_vars( &megno);
	initialize_ResonanceData(&r);
	initialize_ResonanceData(&r1);

	double alphaIn0=pow(1./n1,1./1.5);
	double alphaOut0=pow(n2,1./1.5);

#if INTEGRATE
	int innerRes[3*3];
	int outerRes[3*3];
	for(int j=0;j<3;j++){
		innerRes[3*j]=j+3;
		innerRes[3*j+1]=1;
		innerRes[3*j+2]=0;
		outerRes[3*j]=j+3;
		outerRes[3*j+1]=1;
		outerRes[3*j+2]=1;

	}
	ActionAngleSimulation sim;
	InitializeActionAngleSimulation(&sim, 3,innerRes,3,outerRes,mu1,mu2,n1,n2,e1,e2,varpi2,L0,l0,X0,Y0,dt);
	Nstep = (int)(tFin/dt);

// printf("%.8f\t%.8f\n",sim.parameters.mu1,sim.parameters.mu2);
// printf("%d\t%d\n",sim.rIn.Nres,sim.rOut.Nres);
// printf("%d\t%d\t%d\n",*(sim.rIn.ResonanceIndices),*(sim.rIn.ResonanceIndices+1),*(sim.rIn.ResonanceIndices+2));

	for(int i=0; i<Nstep; i++){	
		if(i%1000==0){
			printf("%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",sim.t,sim.state.L,sim.state.l,sim.state.X,sim.state.Y,sim.megno.megno);
		}
		SimulationStep(&sim);
	}
// 		ActionAngle_H0_Advance( &state , &pars, 0.5 * dt);
// 		ActionAngle_H1_Advance_StormerVerlet(&state, &r , &r1 ,&pars, i*dt ,dt);
// 		ActionAngle_H0_Advance( &state , &pars, 0.5 * dt);
// 		ActionAngle_Get_Var_Dot(&state, &r , &r1 ,&pars, i*dt+dt);
// 		ActionAngle_update_megno_eqations(&state , &megno ,i*dt+dt,dt);

#endif
#if DERIVS
	AddResonance(&r1,7,2,2,alphaOut0 );
	double* derivs =  malloc(4*sizeof(double));
	double* jacobian =  malloc(4*4*sizeof(double));
	H1_Outer_Derivs(derivs,jacobian,&state,&r1,mu1,n1,n2,e2,0,0);

	printf("\nDerivatives\n");
	for(int i=0; i<4; i++){
	printf("%.8g\t",derivs[i]);
	}
	
	printf("\n\nJacobian \n");
	
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			printf("%.5g\t", *(jacobian + INDX(i,j)));
		}
	printf("\n");
	}

	free(derivs);
	free(jacobian);

#endif

	free_ResonanceData(&r);
	free_ResonanceData(&r1);

	return 0;
}
// Check derivative function
