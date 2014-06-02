#include "FR.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int FR::derivatives(const double data[], double derivs[], Parameters &params) {

	// The first section here is model independent

	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];
	double hubblesq = pow(hubble,2.0);
	
	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;
	// Computing \dot{\phi}
	derivs[1] = phidot;

	// Model dependent code below

	// Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
	// (1) Compute FR Lagrangian quantities
	int result = computelagrangian(data);
	
	
	// Split up the acceleration equation: H1 is all the "standard" contributions
	double H1 = - 0.5 * pow(hubble, 2.0) - 0.5 * params.OmegaR() / a2 + 0.5 * params.OmegaK();
	// put it all together and compute \dot{H}
	derivs[3]=H1-1.5*a2*pressure

	// Now that we've computed hubbledot, can finish off the scalar field equation of motion:
	derivs[2] = -1/(2*Fpp)*(a2*P+a2*(F+4*Fppp*X-phi*Fp)+(1+2*Fp)*(hubblesq+2*phidot+curv)+2*Fpp*hubble*phidot);
	
	// To compute the pressure from the function, need to pass it \dot{H}
	double press = pressure(data, derivs[3]);
	
	// GSL_SUCCESS indicates that the computation was successful.
	// If it failed, look up the appropriate error code in the same enum as GSL_SUCCESS
	return GSL_SUCCESS;
}

// Function to calculate the Lagrangian and all its appropriate derivatives
// The results array should be of length 6
int FR::computelagrangian(const double data[]) {
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
	
	// Compute powers of a:
	a2 = pow(data[0],2.0);
	a4 = a2 * a2;
	a8 = a4 * a4;

    F = 1.0;
	Fp = 1.0;
	Fpp = 1.0;
	Fppp = 1.0; 
	
	 
	
	
	 
	
	// Store the data for which these results are correct
	for (int i = 0; i < 4; i++)
		storeddata[i] = data[i];

	// Success!
	return 0;
}


// Returns the ratio rho_Q/rho_c
double FR::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];
	double a2 = pow(a,2.0);
	
	// Compute quantities
	int result = computelagrangian(data);

	// Now construct energy/3 for FR
	return (-(a2*(F-phi*Fp)+6.0*(Fp*pow(hubble,2.0)*Fp*curv+Fpp*hubble*phidot))/a2)/3.0;
}

// Returns the ratio P_Q/rho_c
double FR::pressure(const double data[], const double hdot){
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
	
	// Now put the above bit in for \ddot{\phi}.
	return ()/3.0;
}

/* This function does a few things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
std::string FR::init(double data[], double time, Parameters &params, IniReader &init, int &errorstate) {

	// Set the name of the class
	section = "FR";

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
	
	/// NEED to modify this for computing the initial Hubble for FR
	temp = params.OmegaM() / a + params.OmegaR() / a2 + params.OmegaK() + a2 * energydensity(data);

	// Calculate H
	data[3] = pow(temp, 0.5);

	// We have success!
	errorstate = 0;

	// Return a string to print to the log
	std::stringstream output;
	output << "Running FR model. lambda = " << lambda
			<< ", n = " << n << ", alpha = " << alpha << ", beta = " << beta << std::endl;
	return output.str();

}

// The speedofsound2 returns the speed of sound squared, given the state of the system
double FR::speedofsound2(const double data[]) {
	// The speed of sound in FR is 1.
	// Extract quantities for easier reading of the code
	
		
	// Return result
	return 1.0;
}

// The isghost function is given the state of the system and returns whether or not the theory has become ghostlike
bool FR::isghost(const double data[]) {
	
	// FR becomes ghost-like when the denominator in the speed of sound squared becomes negative.
	double phidot = data[2];
	double hubble = data[3];
	
	// Compute quantities
	int result = computelagrangian(data);

	// Check if positive, and return the result
	if (1+2.0*Fp < 0)
		return true;
	else
		return false;

}
