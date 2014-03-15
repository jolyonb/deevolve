#include "constw.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int ConstW::derivatives(const double data[], double derivs[], Parameters &params) {

	// The first section here is model independent

	// Extract data for easier reading of the code
	double a = data[0];
	double a2 = pow(a, 2.0);   // a^2
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;
	// Computing \dot{\phi}
	derivs[1] = phidot;

	// Model dependent code below

	// Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
	// There is no scalar field in this model, so phidot = 0.
	derivs[2] = 0;

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Note that pressure does not depend on \dot{H} in this model,
	// so we pass in 0 for \dot{H} when calculating pressure
	double press = pressure(data, 0.0);
	derivs[3] = - pow(hubble, 2.0) / 2 - params.OmegaR() / a2 / 2 + params.OmegaK() / 2 - 3.0 * a2 * press / 2.0;

	// GSL_SUCCESS indicates that the computation was successful. If it failed, look up the appropriate error code in the same enum as GSL_SUCCESS
	return GSL_SUCCESS;
}

// Returns the ratio rho_Q/rho_c
double ConstW::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	// Your code here
	return pow(a, - 3.0 * (1.0 + EOSw)) * OmegaLambda;
}
// Returns the ratio P_Q/rho_c
double ConstW::pressure(const double data[], const double hdot){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	// Your code here
	return EOSw * energydensity(data);
}

/* This function does four things:
 * - Sets the name of the class
 * - Reads in the value of Omega_Lambda from the ini file
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
std::string ConstW::init(double data[], double time, Parameters &params, IniReader &init, int &errorstate) {

	// Set the name of the class
	section = "ConstW";

	// Go and extract Omega_Lambda from the ini file
	OmegaLambda = init.getiniDouble("OmegaLambda", 0.7, section);
	// Also extract w from the ini file
	EOSw = init.getiniDouble("EOSw", -1.0, section);

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

	// Discard phi, \dot{\phi}
	data[1] = 0;
	data[2] = 0;

	// We have success!
	errorstate = 0;

	// Return a string to print to the log
	std::stringstream output;
	output << "Running ConstW model with Omega_Lambda = " << OmegaLambda << " and w = " << EOSw << std::endl;
	return output.str();

}
