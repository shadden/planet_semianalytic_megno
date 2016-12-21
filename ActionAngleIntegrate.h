void InitializeActionAngleSimulation(ActionAngleSimulation* sim, int NresIn, int includeInnerZeroth, int* innerRes,int NresOut,int includeOuterZeroth, int* outerRes,\
									 double mu1,double mu2,double n1,double n2,double e1,double e2,double lambda2,double varpi2,\
									 double L0,double l0,double X0, double Y0,double dt);
									 
void SimulationStep(ActionAngleSimulation* restrict sim);
void ActionAnglePhaseStateInitialize(ActionAnglePhaseState* restrict Z, double L0, double l0, double X0, double Y0);
void ActionAngle_H0_Advance( ActionAnglePhaseState* restrict Z ,SimulationParameters* restrict pars,const double t,const double dt);
void initialize_pars(SimulationParameters* pars,double ,double ,double ,double ,double,double,double,double);
void ActionAngle_H1_Advance(ActionAnglePhaseState* Z, ResonanceData* restrict rIn, ResonanceData* restrict rOut ,const double mu1, const double mu2,\
	const double n1, const double n2, const double t, const double dt);
void ActionAngle_H1_Advance_StormerVerlet(ActionAnglePhaseState* Z, ResonanceData* restrict rIn, ResonanceData* restrict rOut ,SimulationParameters* restrict pars, const double t, const double dt);
void ActionAngle_update_megno_eqations(ActionAnglePhaseState* restrict particle , MEGNO_Auxilary_Variables* megno ,double t,double dt);
void ActionAngle_Get_Var_Dot(ActionAnglePhaseState* state, ResonanceData* restrict rIn, ResonanceData* restrict rOut ,SimulationParameters* restrict pars, const double t);

void H1_Inner_Derivs(double* derivs, double* jacobian, ActionAnglePhaseState* Z, ResonanceData* restrict rIn ,const double mu1, const double n1, const double e1, const double t);
void H1_Outer_Derivs(double* derivs, double* jacobian, ActionAnglePhaseState* Z, ResonanceData* restrict rOut ,const double mu2, const double n1,const double n2, const double e2,const double lambda2,const double varpi2, const double t);
