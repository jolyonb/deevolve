#include "kessence.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int Kessence::derivatives(const double data[], double derivs[], Parameters &params) {

	// The first section here is model independent

	// Extract data for easier reading of the code
	double a = data[0];
	double a2 = pow(a, 2.0);   // a^2
	double phi = data[1];
	double phidot = data[2];
	double phidot2 = pow(data[2], 2.0);
	double phidot3 = pow(data[2], 3.0);
	double hubble = data[3];

	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;
	// Computing \dot{\phi}
	derivs[1] = phidot;

	// Model dependent code below

	// Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
	// Compute k-essence quantities
	double U[6];
	int result = computefunctions(data, U);
	double X = pow(data[2] / data[0], 2.0) / 2;

	derivs[2] = (a2 * U[1]
	                   - 2 * hubble * phidot * U[3]
	                   - U[5] * phidot2
	                   + U[4] * hubble * phidot3 / a2)
	                   / (U[3] + U[4] * phidot2 / a2);

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Note that pressure does not depend on \dot{H} in this model,
	// so we pass in 0 for \dot{H} when calculating pressure
	double press = pressure(data, 0.0);
	derivs[3] = - 0.5 * pow(hubble, 2.0) - 0.5 * params.OmegaR() / a2 + 0.5 * params.OmegaK() - 1.5 * a2 * press;

	// GSL_SUCCESS indicates that the computation was successful.
	// If it failed, look up the appropriate error code in the same enum as GSL_SUCCESS
	return GSL_SUCCESS;
}

// Function to calculate the Lagrangian and all its appropriate derivatives:
// U, Up, Upp, UX, UXX, UXP
// The results array should be of length 6
int Kessence::computefunctions(const double data[], double results[]) {
	// First, check the data against the previous data
	// This prevents the results being computed multiple times on the same data
	if (data[0] == storeddata[0] &&
			data[1] == storeddata[1] &&
			data[2] == storeddata[2] &&
			data[3] == storeddata[3]) {
		// It's the same as before, so return the stored results
		for (int i = 0; i < 6; i++)
			results[i] = storedresults[i];
		// Success!
		return 0;
	}

	// Extract data for easier reading of the code
	double phi = data[1];
	double X = pow(data[2] / data[0], 2.0) / 2;

	// Calculate everything
	// U is a function of X and phi
	// Quintessence is U = X - V(phi)
	// This implements quintessence, in order to compare with the quintessence module
	// U
	results[0] = X - pow(phi, 2.0) / 2.0;
	// U,phi
	results[1] = - phi;
	// U,phi phi
	results[2] = -1.0;
	// U,X
	results[3] = 1.0;
	// U,X X
	results[4] = 0.0;
	// U,X phi
	results[5] = 0.0;

	// Store the data and calculated results
	for (int i = 0; i < 4; i++)
		storeddata[i] = data[i];
	for (int i = 0; i < 6; i++)
		storedresults[i] = results[i];

	// Success!
	return 0;
}

// Returns the ratio rho_Q/rho_c
double Kessence::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	// Compute quantities
	double U[6];
	int result = computefunctions(data, U);
	double X = pow(data[2] / data[0], 2.0) / 2;

	return (2 * X * U[3] - U[0]) / 3.0;
}
// Returns the ratio P_Q/rho_c
double Kessence::pressure(const double data[], const double hdot){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	// Compute quantities
	double U[6];
	int result = computefunctions(data, U);

	return U[0] / 3.0;
}

/* This function does a few things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
std::string Kessence::init(double data[], double time, Parameters &params, IniReader &init, int &errorstate) {

	// Set the name of the class
	section = "Kessence";

	// Go and get model parameters
	// lambda is just a parameter in the ini file. It can be used to do a single parameter scan, for example.
	lambda = init.getiniDouble("lambda", 1.0, section);

	// Construct H
	// Temporary variable
	double temp;
	// Scale factor
	double a = data[0];
	double a2 = pow(a, 2.0);

	// Calculate H^2
	temp = params.OmegaM() / a + params.OmegaR() / a2 + params.OmegaK() + a2 * energydensity(data);

	// Calculate H
	data[3] = pow(temp, 0.5);

	// We have success!
	errorstate = 0;

	// Return a string to print to the log
	std::stringstream output;
	output << "Running Kessence model with parameter lambda = " << lambda << std::endl;
	return output.str();

}
