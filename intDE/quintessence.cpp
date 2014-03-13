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


// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi})
int Quintessence::derivatives(const double data[], double derivs[], Parameters &params) {

	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	// Calculate a^2
	double a2 = pow(a, 2.0);
	// Temporary variable
	double temp;
	// Hubble parameter
	double hubble;

	// Computing \dot{a}
	// This one comes from the Friedmann equation
	// temp is \dot{a}^2
	temp = a * params.OmegaM() + params.OmegaR() + a2 * params.OmegaK() + pow(a2, 2.0) * energydensity(data);
	derivs[0] = pow(temp, 0.5);
	hubble = derivs[0] / a;

	// Computing \dot{\phi}
	// This one is easy: \dot{\phi} = \dot{\phi}
	derivs[1] = data[2];

	// Computing \ddot{\phi}
	derivs[2] = - 2.0 * hubble * phidot - a2 * potentialprime(phi);

	// Computing \dot{H}. We don't need that in this model, so we set it to zero.
	derivs[3] = 0;

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
void Quintessence::getstate(const double data[], double time, double info[], Parameters &params) {
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	// Calculate a^2, a^4
	double a2 = pow(a, 2.0);
	double a4 = pow(a, 4.0);
	// Go and compute the derivatives
	double derivs[4];
	int result = derivatives(data, derivs, params);
	// Hubble parameter H = \dot{a}/a
	double hubble = derivs[0] / a;
	double hubble2 = pow(hubble, 2.0);
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

	// phi
	info[5] = phi;

	// \dot{phi}
	info[6] = phidot;

	// \ddot{phi}
	info[7] = derivs[2];

	// rho_Q / rho_c
	info[13] = energy;

	// P_Q / rho_c
	info[14] = press;

	// w_Q
	info[15] = press / energy;

	// \dot{H}
	info[4] = - params.OmegaM() / 2 / a - params.OmegaR() / a2 - a2 * (3.0 * press + energy) / 2.0;

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

	// Error
	info[16] = 0;

}

// Returns the ratio rho_Q/rho_c
double Quintessence::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];

	return (pow(phidot,2.0) / 2 / pow(a,2.0) + potential(phi)) / 3;
}
// Returns the ratio P_Q/rho_c
double Quintessence::pressure(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];

	return (pow(phidot,2.0) / 2 / pow(a,2.0) - potential(phi)) / 3;
}

// The speedofsound2 returns the speed of sound squared, given the state of the system
double Quintessence::speedofsound2(const double data[]) {
	// The speed of sound in quintessence is always 1.
	return 1;
}
// The implementsSOS function returns whether or not a class actually implements the speedofsound2 function
bool Quintessence::implementsSOS() {
	return true;
}

// The isghost function is given the state of the system (and some derivatives, because those might be helpful)
// and returns whether or not the theory has become ghostlike
bool Quintessence::isghost(const double data[]) {
	// Quintessence is never ghost-like
	return false;
}
// The implementsghost function returns whether or not a class actually implements the isghost function
bool Quintessence::implementsghost() {
	return true;
}

// Return the description of the model
std::string Quintessence::description() {
	std::stringstream output;
	output << "Running Quintessence model." << std::endl;
	return output.str();
}
