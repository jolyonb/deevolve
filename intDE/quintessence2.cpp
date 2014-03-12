#include "quintessence2.h"

// To modify this module to a given quintessence potential, the following two functions must be specified
// Returns the potential for a given phi (used in calculating energy density and pressure)
double QuintessenceH::potential(const double phi){

	// This is just the potential V = m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
	// In terms of the dimensionless potential U, this is U = phi^2 / 2
	return pow(phi, 2.0) / 2;

	// This is an exponential function: V = m_P^2 H_0^2 e^(-lambda phi)
	//const double lambda = 1.0;
	//return exp(- lambda * phi);
}
// Returns the derivative of the potential for a given phi (used in calculating the scalar equation of motion)
double QuintessenceH::potentialprime(const double phi){

	// This is just the potential V = m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
	// In terms of the dimensionless potential U, this is U = phi^2 / 2
	// \partial_\phi U = phi, so return phi
	return phi;

	// This is an exponential function: V = m_P^2 H_0^2 e^(-lambda phi)
	//const double lambda = 1.0;
	//return - lambda * exp(- lambda * phi);

}

// This function initializes the value of H, using the Friedmann equation
int QuintessenceH::init(double data[], double time, Parameters &params) {

	// Temporary variable
	double temp;
	// Scale factor
	double a = data[0];
	double a2 = pow(a, 2.0);

	// Calculate H^2
	temp = params.OmegaM() / a + params.OmegaR() / a2 + params.OmegaK() + a2 * energydensity(data);

	// Calculate H
	data[3] = pow(temp, 0.5);

	return 0; // Success!
}

// The derivatives routine is given the state of the system (a, phi and \dot{\phi} as well as H in some cases)
// as well as the parameters of the system, and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi} (and \dot{H}, if necessary))
// In this module, we are evolving a by integrating H.
int QuintessenceH::derivatives(const double data[], double derivs[], Parameters &params) {

	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];
	// Calculate a^2
	double a2 = pow(a, 2.0);
	// Expressions for energy and pressure of dark energy
	double press = pressure(data);
	double energy = energydensity(data);

	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;

	// Computing \dot{\phi}
	// This one is easy: \dot{\phi} = \dot{\phi}
	derivs[1] = phidot;

	// Computing \ddot{\phi}
	derivs[2] = - 2.0 * hubble * phidot - a2 * potentialprime(phi);

	// Computing \dot{H}. This requires the acceleration equation.
	derivs[3] = - params.OmegaM() / 2 / a - params.OmegaR() / a2 - a2 * (3.0 * press + energy) / 2.0;

	return GSL_SUCCESS;
}

/* The state routine is given the state of the system as well as the parameters of the model,
   and returns information in the info array. The return values are as follows:

   * 0 time
   * 1 a
   * 2 Redshift
   * 3 H = \dot{a}/a
   * 4 \dot{H}
   * 5 phi
   * 6 \dot{phi}
   * 7 \ddot{phi}
   * 8 Omega_matter (present value)
   * 9 Omega_radiation (present value)
   * 10 Omega_k (present value)
   * 11 Omega_Q (present value)
   * 12 w_total
   * 13 rho_Q / rho_c
   * 14 P_Q / rho_c
   * 15 w_Q
   * 16 Error

 */
void QuintessenceH::getstate(const double data[], double time, double info[], Parameters &params) {
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];
	double hubble2 = pow(hubble, 2.0);
	// Calculate a^2, a^4
	double a2 = pow(a, 2.0);
	double a4 = pow(a, 4.0);
	// Go and compute the derivatives
	double derivs[4];
	int result = derivatives(data, derivs, params);
	// Energy density and pressure
	double energy = energydensity(data);
	double press = pressure(data);

	// time
	info[0] = time;

	// a
	info[1] = a;

	// Redshift
	info[2] = 1.0 / a - 1.0;

	// Hubble
	info[3] = hubble;

	// \dot{H}
	info[4] = derivs[3];

	// phi
	info[5] = phi;

	// \dot{phi}
	info[6] = phidot;

	// \ddot{phi}
	info[7] = derivs[2];

	// Omega_matter (present value)
	info[8] = params.OmegaM() / a / hubble2;

	// Omega_radiation (present value)
	info[9] = params.OmegaR() / a2 / hubble2;

	// Omega_k (present value)
	info[10] = params.OmegaK() / hubble2;

	// Omega_Q (present value)
	info[11] = a2 / hubble2 * energy;

	// w_total
	info[12] = (params.OmegaR() / 3 + a4 * press) / (a * params.OmegaM() + params.OmegaR() + a4 * energy);

	// rho_Q / rho_c
	info[13] = energy;

	// P_Q / rho_c
	info[14] = press;

	// w_Q
	info[15] = press / energy;

	// Error
	// What is the relative error in our Hubble constant? This expression substitutes the hubble parameter into
	// the Friedmann equation with everything on one side; the other side should thus be zero.
	// The error is expressed as a relative error.
	info[16] = (hubble - pow(params.OmegaM() / a + params.OmegaR() / a2 + params.OmegaK() + a2 * energydensity(data), 0.5) ) / hubble;

}

// Returns the ratio rho_Q/rho_c
double QuintessenceH::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];

	return (pow(phidot,2.0) / 2 / pow(a,2.0) + potential(phi)) / 3;
}
// Returns the ratio P_Q/rho_c
double QuintessenceH::pressure(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];

	return (pow(phidot,2.0) / 2 / pow(a,2.0) - potential(phi)) / 3;
}

// The speedofsound2 returns the speed of sound squared, given the state of the system
double QuintessenceH::speedofsound2(const double data[], const double derivs[]) {
	// The speed of sound in quintessence is always 1.
	return 1;
}
// The implementsSOS function returns whether or not a class actually implements the speedofsound2 function
bool QuintessenceH::implementsSOS() {
	return true;
}

// The isghost function is given the state of the system (and some derivatives, because those might be helpful)
// and returns whether or not the theory has become ghostlike
bool QuintessenceH::isghost(const double data[], const double derivs[]) {
	// Quintessence is never ghost-like
	return false;
}
// The implementsghost function returns whether or not a class acutlaly implements the isghost function
bool QuintessenceH::implementsghost() {
	return true;
}
