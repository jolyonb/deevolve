#include "k-essence.h"

// To modify this module to a given quintessence potential, the following two functions must be specified
// Returns the potential for a given phi (used in calculating energy density and pressure)
/*
double Kessence::potential(const double phi){

	// This is just the potential V = m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
	// In terms of the dimensionless potential U, this is U = phi^2 / 2
	return pow(phi, 2.0) / 2;

	// This is an exponential function: V = m_P^2 H_0^2 e^(-lambda phi)
	//const double lambda = 1.0;
	//return exp(- lambda * phi);
}
// Returns the derivative of the potential for a given phi (used in calculating the scalar equation of motion)
double Kessence::potentialprime(const double phi){

	// This is just the potential V = m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
	// In terms of the dimensionless potential U, this is U = phi^2 / 2
	// \partial_\phi U = phi, so return phi
	return phi;

	// This is an exponential function: V = m_P^2 H_0^2 e^(-lambda phi)
	//const double lambda = 1.0;
	//return - lambda * exp(- lambda * phi);

}
*/
double* Kessence::LagrangianQuants(const double data[]){
    
    double a = data[0];
    double phi = data[1];
    double phidot = data[2];
    double a2 = pow(a,2.0);
    
    // Array which holds info about the Lagrangian
    double *LagQuants[6];
    
// BEGIN: specific choice of parameterization of the Lagrangian
    // This Lagrangian is L = X^n - V(phi)
    //  with V = e^{-lambda * phi}
    //   (all suitably dimensionless)
    
    // Power in exponential in potential
    const double lambda = 2.0;
    
    // Power of kinetic term
    // When n = 1, should reproduce quintessence
    const double n = 1.0;
    
    // Here is the potential
    double pot = exp(- lambda * phi);
    double dpot = - lambda * pot;
    
    // Here we construct the Lagrangian & important derivatives thereof.
    
    // (1) Kinetic scalar
    double X = pow( phidot , 2.0 ) / a2 / 2.0;
    
    // (2) Lagrangian
    double L = pow( X, n ) - pot;
    
    // (3) dL/dX
    double LX = n * pow ( X, n - 1.0 );
    
    // (4) d^2L/dX^2
    double LXX = n * (n - 1.0) * pow( X, n - 2.0);
    
    // (5) dL/dphi
    double Lp = - dpot;
    
    // (6) d^2L/dXdphi
    double LpX = 0.0;

// END: specific choice of parameterization of the Lagrangian
    
    // Put these choices into LagQuants array, to be passed back
    
    LagQuants[0] = X;
    LagQuants[1] = L;
    LagQuants[2] = LX;
    LagQuants[3] = LXX;
    LagQuants[4] = Lp;
    LagQuants[5] = LpX;
    
    return LagQuants;
    
}

// This function initializes the value of H, using the Friedmann equation
int Kessence::init(double data[], double time, Parameters &params) {

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
int Kessence::derivatives(const double data[], double derivs[], Parameters &params) {

	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];
	// Calculate a^2
	double a2 = pow(a, 2.0);
    
    // Quantities about the Lagrangian (derivatives, potential etc)
    double *LagQuants = LagrangianQuants(data);
    
	// Expressions for energy, pressure, sound speed of dark energy
	double press = pressure(LagQuants);
	double energy = energydensity(LagQuants);
    double SoS = speedofsound2(LagQuants);
    
    double X = LagQuants[0];
    double L = LagQuants[1];
    double LX = LagQuants[2];
//  double LXX = LagQuants[3];      // not actually used in the EoM (inside cs2)
    double Lp = LagQuants[4];
    double LpX = LagQuants[5];

    
	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;

	// Computing \dot{\phi}
	// This one is easy: \dot{\phi} = \dot{\phi}
	derivs[1] = phidot;

	// Computing \ddot{\phi}
	derivs[2] = ( 1.0 - 3.0 * SoS ) * phidot * hubble + SoS / LX * ( a2 * Lp - pow( phidot , 2.0 ) * LpX );

	// Computing \dot{H}. This requires the acceleration equation.
	derivs[3] = - 0.5 * pow( hubble , 2.0 ) - params.OmegaR() / 2.0 / a2 - 1.5 * a2 *  L;

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

void Kessence::getstate(const double data[], double time, double info[], Parameters &params) {
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
    
    double *LagQuants = LagrangianQuants(data);

    // Energy density, pressure, and SoS
	double energy = energydensity(LagQuants);
	double press = pressure(LagQuants);
    double SoS = speedofsound2(LagQuants);
    
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
double Kessence::energydensity(const double LagQuants[]){
	// Extract data for easier reading of the code
    
	double X = LagQuants[0];
	double L = LagQuants[1];
	double LX = LagQuants[2];

	return ( 2.0 * X * LX - L ) / 3.0;
}
// Returns the ratio P_Q/rho_c
double Kessence::pressure(const double LagQuants[]){
	// Extract data for easier reading of the code
    
	double L = LagQuants[1];
    
	return L / 3.0;
}

// The speedofsound2 returns the speed of sound squared, given the state of the system
double Kessence::speedofsound2(const double LagQuants[]) {
    
    double X = LagQuants[0];
    double LX = LagQuants[2];
    double LXX = LagQuants[3];
    
	return LX / ( LX + 2.0 * X * LXX );

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
