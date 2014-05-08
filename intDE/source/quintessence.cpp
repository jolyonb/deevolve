#include "quintessence.h"

// To modify this module to a given quintessence potential, the following two functions must be specified
// Returns the potential for a given phi (used in calculating energy density and pressure)
double Quintessence::potential(const double phi){

	// The potential is selected by a parameter that is specified in params.ini
	switch(modeltype){
		case 1 :
			// lambda phi^4
			// This is just the potential V = mass * m_P^2 H_0^2 phi^2 / 2 + lambda * m_P^2 H_0^2 phi^4
			// (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = mass * phi^2 / 2 + lambda * phi^4
			return mass * pow(phi, 2.0) / 2 + lambda * pow(phi, 4.0);

		case 2 :
			// Exponential
			// This is an exponential function: V = alpha m_P^2 H_0^2 exp(- beta phi)
			// (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = alpha exp(- beta phi)
			return alpha * exp(- beta * phi);

		case 3 :
			// User defined potential
			// (Presently not defined)
			return 0;

		case 0 :
		default : // Catch all
			// mass term only
			// This is just the potential V = mass * m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = mass * phi^2 / 2
			return mass * pow(phi, 2.0) / 2;
	}

	// It's an error if we get to here, but return something sensible
	return 0;

}
// Returns the derivative of the potential for a given phi (used in calculating the scalar equation of motion)
double Quintessence::potentialprime(const double phi){

	// The potential is selected by a parameter that is specified in params.ini
	switch(modeltype){
		case 1 :
			// lambda phi^4
			// This is just the potential V = mass * m_P^2 H_0^2 phi^2 / 2 + lambda * m_P^2 H_0^2 phi^4
			// (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = mass * phi^2 / 2 + lambda * phi^4
			// The derivative is U' = mass * phi + 4 * lambda * phi^3
			return mass * phi + 4.0 * lambda * pow(phi, 3.0);

		case 2 :
			// Exponential
			// This is an exponential function: V = alpha m_P^2 H_0^2 exp(- beta phi)
			// (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = alpha exp(- beta phi)
			// The derivative is U' = - alpha beta exp(- beta phi)
			return - beta * alpha * exp(- beta * phi);

		case 3 :
			// User defined potential
			// (Presently not defined)
			return 0;

		case 0 :
		default : // Catch all
			// mass term only
			// This is just the potential V = mass * m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = mass * phi^2 / 2
			// The derivative is U' = mass * phi
			return mass * phi;
	}

	// It's an error if we get to here, but return something sensible
	return 0;

}

/* This function does three things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
std::string Quintessence::init(double data[], double time, Parameters &params, IniReader &init, int &errorstate) {

	// Set class name
	section = "Quintessence";

	// Go and get model parameters
	modeltype = init.getiniInt("PotentialType", 0, section);
	if (modeltype < 0 || modeltype > 3)  // do a quick bit of error checking
		modeltype = 0;
	mass = init.getiniDouble("mass", 1.0, section);
	lambda = init.getiniDouble("lambda", 1.0, section);
	alpha = init.getiniDouble("alpha", 1.0, section);
	beta = init.getiniDouble("beta", 1.0, section);

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

	// Mark success!
	errorstate = 0;

	// Return the description of the model
	std::stringstream output;
	switch(modeltype){
		case 1 :
			// lambda phi^4
			output << "Running Quintessence model with lambda phi^4 potential. Mass = "
					<< mass << ", lambda = " << lambda << std::endl;
			break;
		case 2 :
			// Exponential
			output << "Running Quintessence model with exponential potential. alpha = "
					<< alpha << ", beta = " << beta << std::endl;
			break;
		case 3 :
			// User defined potential
			output << "Running Quintessence model with user-defined potential." << std::endl;
			break;
		case 0 :
		default : // Catch all
			// mass term only
			output << "Running Quintessence model with mass potential. Mass = "
					<< mass << std::endl;
	}
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
