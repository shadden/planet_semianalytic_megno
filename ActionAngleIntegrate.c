// Evolve under action of Hamiltonian
// H0 = 0.5 (X^2 + Y^2) - L - 4/L^2
// i.e. Kelerian Hamiltonian after transforming angle
// va
void ActionAngle_H0_Advance( ActionAnglePhaseState* Z , double dt){
		
	const ActionAnglePhaseState state = *Z;
	// Equations
	(*Z.l) += dt * ( 8./(state.L*state.L*state.L) - 1. );
	(*Z.X)  = cos(dt) * state.X - sin(dt) * state.Y;
	(*Z.Y)  = sin(dt) * state.X + cos(dt) * state.Y;
	// Variationals
	(*Z.dl) += dt *  -24./(state.L*state.L*state.L*state.L) * state.dL ;
	(*Z.dX)  = cos(dt) * state.dX - sin(dt) * state.dY;
	(*Z.dY)  = sin(dt) * state.dX + cos(dt) * state.dY;
}

ActionAngle_H1_Advance()
