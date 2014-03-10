#include "params.h"

// The constructor sets the parameters of the model, which cannot be changed afterwards
Parameters::Parameters(const double OmegaM, const double Tgamma, const double OmegaK, const double z0, const double h) {

	const double hubble0conv = 3.24e-18; // This is 100 km/s/Mpc in units of s^-1

	mOmegaM = OmegaM;
	mOmegaK = OmegaK;
	mT = Tgamma;
	mh = h;
	mz0 = z0;
	mH0 = hubble0conv * h;

	// TO DO: This is just a random value until I can find the appropriate formula
	mOmegaR = 0.00001;
	mrhoc = 0;

}
