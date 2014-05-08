#include "KGB.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int KGB::derivatives(const double data[], double derivs[], Parameters &params) {

	// The first section here is model independent

	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double phidot2 = phidot * phidot;
	double phidot3 = phidot * phidot2;
	double phidot4 = phidot2 * phidot2;
	double hubble = data[3];

	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;
	// Computing \dot{\phi}
	derivs[1] = phidot;

	// Model dependent code below

	// Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
	// (1) Compute KGB Lagrangian quantities
	int result = computelagrangian(data);
	
	// Compute \ddot{\phi} = N/D,
	// and split the numerator N and denominator D
	// into contributions from V & L3 separately.
	// N = NV + NL3, D = DV + DL3.
	
	// Numerator for the "V"-part
	double NV = a8*(a2*Vp+hubble*(VXX*phidot3/a2-2.0*VX*phidot)-VXp*phidot2); 
	// Numerator for L3 part; since NL3 contains hubbledot, separate it out 
	// NL3 = NL3_1 + NL3_2 * hubble
	double NL3_1 = a4*(4.0*a4*hubble*phidot*L3p
		+3.0*pow(hubble,2.0)*phidot4*L3XX*a4*phidot2*L3pp-4.0*a2*hubble*phidot3*L3X);
	double NL3_2 = -a4*3.0*a2*phidot2*L3X;
	// Denominator for "V" part
	double DV = a8*(VXX*phidot2/a2 + VX);
	// denominator for "L3" part
	double DL3 = a4*(-2.0*a4*L3p+6.0*a2*hubble*phidot*L3X-a2*phidot2*L3X+3.0*L3XX*hubble*phidot3);
	// Full denominator
	double D = DV + DL3;
	
	// Finally, construct \ddot{phi} = N / D; for that we need hubbledot,
	// and so we wait.

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Note that pressure does depend on \dot{H} in this model,
	
	// Write the pressure as P = P1 + P2\ddot(phi)
	double P1 = (V-2.0*X/a2*(a2*L3p-L3X*hubble*phidot))/3.0;
	double P2 = -2.0*X/a2*L3X/3.0;
	// pressure = P1 + P2\ddot{phi}.

	// Split up the acceleration equation: H1 is all the standard contributions
	double H1 = - 0.5 * pow(hubble, 2.0) - 0.5 * params.OmegaR() / a2 + 0.5 * params.OmegaK();
	// put it all together and compute \dot{H}
	derivs[3]=pow(1.0+1.5*a2*P2*NL3_2/D,-1.0) * (H1-1.5*a2*(P1+P2*(NV + NL3_1)/D));

	// Now that we've computed hubbledot, can finish off the scalar field equation of motino:
	derivs[2] = ( NV + NL3_1 + NL3_2 * derivs[3] ) / D;
	
	// To compute the pressure from the function, need to pass it \dot{H}
	double press = pressure(data, derivs[3]);
	
	// GSL_SUCCESS indicates that the computation was successful.
	// If it failed, look up the appropriate error code in the same enum as GSL_SUCCESS
	return GSL_SUCCESS;
}

// Function to calculate the Lagrangian and all its appropriate derivatives:
// U, Up, Upp, UX, UXX, UXP
// The results array should be of length 6
int KGB::computelagrangian(const double data[]) {
	// First, check the data against the previous data
	// This prevents the results being computed multiple times on the same data
	if (data[0] == storeddata[0] &&
			data[1] == storeddata[1] &&
			data[2] == storeddata[2] &&
			data[3] == storeddata[3]) {
		// It's the same as before, so don't recompute it
		// Success!
		return 0;
	}

	// Extract data for easier reading of the code
	double phi = data[1];
	// Compute the X quantity
	X = pow(data[2] / data[0], 2.0) / 2.0;

	// Calculate everything
	// Below is a different choice of Lagrangian
    // This Lagrangian is L = X + alpha X^n - beta e^{-lambda phi}
    // (all suitably dimensionless)
	// Quintessence has alpha = 0

    // Here is the potential
    double pot = beta * exp(- lambda * phi);
    double dpot = - lambda * pot;

	double a2 = pow(data[0],2.0);
	double a4 = a2 * a2;
	double a8 = a4 * a4;

    // Here we construct the Lagrangian & important derivatives thereof.

    // Construct the two functions V(phi, X) & L3(phi, X).
    V = 1.0;
	L3 = 1.0;
	
	// Now construct their derivatives.
	
	VX = 1.0;
	VXX = 1.0;

	L3p = 1.0;
	L3pp = 1.0;
	L3X = 1.0;
	L3XX = 1.0;
	L3Xp = 1.0;
	
//	L = V - L3 * boxphi;


	//  Useful parameterizations:
	// V = c_2X^n,
	// L3 = c_3 X^m
	// .. "generalized galileon"
	// Could also easily add in potential terms.
	
	
	// Store the data for which these results are correct
	for (int i = 0; i < 4; i++)
		storeddata[i] = data[i];

	// Success!
	return 0;
}

// Returns the ratio rho_Q/rho_c
double KGB::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];
	double a2 = pow(a,2.0);
	
	// Compute quantities
	int result = computelagrangian(data);
	
	// for ease, compute rho for the k-essence case
	double rhok=2.0*X*VX-V;
	
	// Now construct energy/3 for KGB
	return (rhok+2.0*X/a2*(-a2*L3p+3.0*L3X*hubble*phidot))/3.0;
}
// Returns the ratio P_Q/rho_c
double KGB::pressure(const double data[], const double hdot){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double phidot2 = phidot * phidot;
	double phidot3 = phidot * phidot2;
	double phidot4 = phidot2 * phidot2;
	double hubble = data[3];
	
	// Compute quantities
	int result = computelagrangian(data);


	// To compute the pressure, we have to use hdot (i.e. \dot{H}) 
	// instead of \ddot{\phi}.
	// First, compute \ddot{phi} = N / D,
	// with N = NV + NL3_1 + NL3_2*hdot,
	// and D = DV + DL3.
	
	// Numerator for the "V"-part
	double NV = a8*(a2*Vp+hubble*(VXX*phidot3/a2-2.0*VX*phidot)-VXp*phidot2); 
	// Numerator for L3 part; since NL3 contains hubbledot, separate it out 
	// NL3 = NL3_1 + NL3_2 * hubble
	double NL3_1 = a4*(4.0*a4*hubble*phidot*L3p
		+3.0*pow(hubble,2.0)*phidot4*L3XX*a4*phidot2*L3pp-4.0*a2*hubble*phidot3*L3X);
	double NL3_2 = -a4*3.0*a2*phidot2*L3X;
	// Denominator for "V" part
	double DV = a8*(VXX*phidot2/a2 + VX);
	// denominator for "L3" part
	double DL3 = a4*(-2.0*a4*L3p+6.0*a2*hubble*phidot*L3X-a2*phidot2*L3X+3.0*L3XX*hubble*phidot3);
	// Full denominator
	double D = DV + DL3;
	// Hence, we have now computed
	// phidotdot = (NV + NL3_1 + NL3_2*hdot) / D.
	
	// Write the pressure as P = P1 + P2\ddot(phi)
	double P1 = (V-2.0*X/a2*(a2*L3p-L3X*hubble*phidot))/3.0;
	double P2 = -2.0*X/a2*L3X/3.0;
	// Now put the above bit in for \ddot{\phi}.
	return P1 + P2 * (NV + NL3_1 + NL3_2 * hdot) / D;
}

/* This function does a few things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
std::string KGB::init(double data[], double time, Parameters &params, IniReader &init, int &errorstate) {

	// Set the name of the class
	section = "KGB";

	// Go and get model parameters
	// lambda is just a parameter in the ini file. It can be used to do a single parameter scan, for example.
	lambda = init.getiniDouble("lambda", 1.0, section);
	alpha = init.getiniDouble("alpha", 1.0, section);
	beta = init.getiniDouble("beta", 1.0, section);
	n = init.getiniDouble("n", 2.0, section);

	// Construct H
	// Temporary variable
	double temp;
	// Scale factor
	double a = data[0];
	double a2 = pow(a, 2.0);

	// Calculate H^2
	
	/// NEED to modify this for computing the initial Hubble for KGB
	temp = params.OmegaM() / a + params.OmegaR() / a2 + params.OmegaK() + a2 * energydensity(data);

	// Calculate H
	data[3] = pow(temp, 0.5);

	// We have success!
	errorstate = 0;

	// Return a string to print to the log
	std::stringstream output;
	output << "Running KGB model. lambda = " << lambda
			<< ", n = " << n << ", alpha = " << alpha << ", beta = " << beta << std::endl;
	return output.str();

}

// The speedofsound2 returns the speed of sound squared, given the state of the system
double KGB::speedofsound2(const double data[]) {
	// The speed of sound in KGB can vary from 1.
	// Extract quantities for easier reading of the code
	
	double phidot = data[2];
	double hubble = data[3];
	
	
	// Compute quantities
	int result = computelagrangian(data);

	// Compute numerator of cs2
	double numer = a2*(-2.0*L3p+VX+2.0*L3Xp*X-2.0*pow(X*L3X,2.0));
	numer+= 2.0*hubble*(L3X-X*L3XX)*phidot;
	
	// NOTE: this cs2 needs \ddot{\phi}
	//numer+=2.0*(L3X+X*L3XX)*phidotdot;
	// Compute denominator of cs2
	double denom=a2*(-2.0*L3p+VX+2.0*X*(VXX-L3Xp)+6.0*pow(X*L3X,2.0));
	denom+= 6.0*hubble*(L3X+X*L3XX) * phidot;

	// Return result
	return numer / denom;
}

// The isghost function is given the state of the system and returns whether or not the theory has become ghostlike
bool KGB::isghost(const double data[]) {
	
	// KGB becomes ghost-like when the denominator in the speed of sound squared becomes negative.
	double phidot = data[2];
	double hubble = data[3];
	
	// Compute quantities
	int result = computelagrangian(data);

	
	// Compute quantity which must be positive for stability
	double stab=a2*(-2.0*L3p+VX+2.0*X*(VXX-L3Xp)+6.0*pow(X*L3X,2.0));
	stab+=6.0*hubble*(L3X+X*L3XX)*phidot;

	// Check if positive, and return the result
	if (stab < 0)
		return true;
	else
		return false;

}
