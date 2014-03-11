/*
 * params.h
 *
 * This is a storage class for cosmological parameters. It is given a variety of inputs, computes some results from that,
 * and then just holds them for future reference.
 *
 */

#ifndef PARAMS_H_
#define PARAMS_H_

#include <cmath>

class Parameters {
	public:
		Parameters(const double OmegaM, const double Tgamma, const double OmegaK, const double z0, const double h); // Constructor

		// Getters for the parameters of the model
		inline double OmegaM () {return mOmegaM;}  // Energy fraction of matter today
		inline double OmegaR () {return mOmegaR;}  // Energy fraction of radiation today
		inline double OmegaK () {return mOmegaK;}  // Energy fraction of spatial curvature today
		inline double Tgamma () {return mT;}       // Temperature of photons today
		inline double h () {return mh;}            // Hubble parameter is h * 100 km/s/Mpc
		inline double H0 () {return mH0;}          // Hubble parameter today (in s^-1)
		inline double z0 () {return mz0;}          // Redshift to start the evolution at
		inline double rhoc () {return mrhoc;}      // Critical energy density (today), in units of eV^4

	private:
		// Internal storage for values
		double mOmegaM;
		double mOmegaR;
		double mOmegaK;
		double mT;
		double mh;
		double mH0;
		double mz0;
		double mrhoc;

};

#endif /* PARAMS_H_ */
