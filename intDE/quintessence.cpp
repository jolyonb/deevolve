#include "quintessence.h"

// To modify this module to a given quintessence potential, the following two functions must be specified
// Returns the potential for a given phi (used in calculating energy density and pressure)
double Quintessence::potential(const double phi){

	// This is just the potential V = m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
	// In terms of the dimensionless potential U, this is U = phi^2 / 2
	return pow(phi, 2.0) / 2;

	// This is an exponential function: V = m_P^2 H_0^2 e^(-lambda phi)
	//const double lambda = 1.0;
	//return exp(- lambda * phi);
}
// Returns the derivative of the potential for a given phi (used in calculating the scalar equation of motion)
double Quintessence::potentialprime(const double phi){

	// This is just the potential V = m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
	// In terms of the dimensionless potential U, this is U = phi^2 / 2
	// \partial_\phi U = phi, so return phi
	return phi;

	// This is an exponential function: V = m_P^2 H_0^2 e^(-lambda phi)
	//const double lambda = 1.0;
	//return - lambda * exp(- lambda * phi);

}

/* This function does three things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
std::string Quintessence::init(double data[], double time, Parameters &params, IniReader &init, int &errorstate) {

	// Set class name
	section = "Quintessence";

	// Temporary variable
	double temp;
	// Scale factor
	double a = data[0];
	double a2 = pow(a, 2.0);

	// Calculate H^2
	temp = params.OmegaM() / a + params.OmegaR() / a2 + params.OmegaK() + a2 * energydensity(data);

	// Calculate H
	data[3] = pow(temp, 0.5);

	// Mark success!
	errorstate = 0;

	// Return the description of the model
	std::stringstream output;
	output << "Running Quintessence model." << std::endl;
	return output.str();

}

/*
 * The derivatives routine is given the state of the system (a, phi, \dot{\phi} and H)
 * as well as the parameters of the system, and returns the derivatives
 * \dot{a}, \dot{\phi}, \ddot{\phi} and \dot{H}
 */
int Quintessence::derivatives(const double data[], double derivs[], Parameters &params) {

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
	derivs[2] = - 2.0 * hubble * phidot - a2 * potentialprime(phi);

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Expressions for energy and pressure of dark energy
	// Note that pressure does not depend on \dot{H} in this model,
	// so we pass in 0 for \dot{H} when calculating pressure
	double press = pressure(data, 0.0);
	derivs[3] = - 0.5 * pow(hubble, 2.0) - 0.5 * params.OmegaR() / a2 + 0.5 * params.OmegaK() - 1.5 * a2 * press;

	return GSL_SUCCESS;
}

// Returns the ratio rho_Q/rho_c
double Quintessence::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	return (0.5 * pow(phidot / a, 2.0) + potential(phi)) / 3;
	// The factor of 1/3 is correct. Note that a cosmological constant will contribute
	// 8 pi G Lambda / 3 H_0^2 = Lambda / \rho_c = Omega_Lambda.
}
// Returns the ratio P_Q/rho_c
double Quintessence::pressure(const double data[], const double hdot){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	return (0.5 * pow(phidot / a, 2.0) - potential(phi)) / 3;
}

// The implementsSOS function returns whether or not a class actually implements the speedofsound2 function
bool Quintessence::implementsSOS() {
	return true;
}
// The speedofsound2 returns the speed of sound squared, given the state of the system
double Quintessence::speedofsound2(const double data[]) {
	// The speed of sound in quintessence is always 1.
	return 1.0;
}

// The implementsghost function returns whether or not a class actually implements the isghost function
bool Quintessence::implementsghost() {
	return true;
}
// The isghost function is given the state of the system and returns whether or not the theory has become ghostlike
bool Quintessence::isghost(const double data[]) {
	// Quintessence is never ghost-like
	return false;
}
