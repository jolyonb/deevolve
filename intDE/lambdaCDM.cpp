#include "lambdaCDM.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi})
int LambdaCDM::derivatives(const double data[], double derivs[], Parameters &params) {

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
	derivs[2] = 0;

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Expressions for energy and pressure of dark energy
	// Note that pressure does not depend on \dot{H} in this model,
	// so we pass in 0 for \dot{H} when calculating pressure
	double press = pressure(data, 0.0);
	// double energy = energydensity(data);
	// derivs[3] = - params.OmegaM() / 2 / a - params.OmegaR() / a2 - a2 * (3.0 * press + energy) / 2.0;
	derivs[3] = - pow(hubble, 2.0) / 2 - params.OmegaR() / a2 / 2 + params.OmegaK() / 2 - 3.0 * a2 * press / 2.0;

	return GSL_SUCCESS;
}

// Returns the ratio rho_Q/rho_c
double LambdaCDM::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	return OmegaLambda;
}
// Returns the ratio P_Q/rho_c
double LambdaCDM::pressure(const double data[], const double hdot){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	return -OmegaLambda;
}

/* This function does four things:
 * - Sets the name of the class
 * - Reads in the value of Omega_Lambda from the ini file
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
std::string LambdaCDM::init(double data[], double time, Parameters &params, IniReader &init, int &errorstate) {

	// Set the name of the class
	section = "LambdaCDM";

	// Go and extract Omega_Lambda = Lambda 8 pi G / 3 / H_0^2 from the ini file
	// where S = \int d^4x \sqrt{-g} ( m_P^2/2 R - Lambda ) defines Lambda
	OmegaLambda = init.getiniDouble("OmegaLambda", 0.7, section);

	// Construct H
	// Temporary variable
	double temp;
	// Scale factor
	double a = data[0];
	double a2 = pow(a, 2.0);

	// Calculate H^2
	temp = params.OmegaM() / a + params.OmegaR() / a2 + params.OmegaK() + a2 * OmegaLambda;

	// Calculate H
	data[3] = pow(temp, 0.5);

	// Discard phi, \dot{\phi}
	data[1] = 0;
	data[2] = 0;

	// We have success!
	errorstate = 0;

	// Return a string to print to the log
	std::stringstream output;
	output << "Running LambdaCDM model with Omega_Lambda = " << OmegaLambda << std::endl;
	return output.str();

}
