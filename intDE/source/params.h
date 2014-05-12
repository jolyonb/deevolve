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
#include "inireader.h"

class Parameters {
	public:

		// Constructor. Sets the parameters of the model, which cannot be changed afterwards
		Parameters (IniReader&);

		// Getters for the parameters of the model
		inline double OmegaM () {return mOmegaM;}  // Energy fraction of matter today
		inline double OmegaB () {return mOmegaB;}  // Energy fraction of baryonic matter today
		inline double OmegaR () {return mOmegaR;}  // Energy fraction of radiation today
		inline double OmegaGamma () {return mOmegaGamma;}  // Energy fraction of photons today
		inline double OmegaK () {return mOmegaK;}  // Energy fraction of spatial curvature today
		inline double Tgamma () {return mT;}       // Temperature of photons today
		inline double h () {return mh;}            // Hubble parameter is h * 100 km/s/Mpc
		inline double H0 () {return mH0;}          // Hubble parameter today (in s^-1)
		inline double Neff () {return mNeff;}      // Effective number of relativistic species (3.046 for standard neutrinos)
		inline double z0 () {return mz0;}          // Redshift to start the evolution at
		inline double rhoc () {return mrhoc;}      // Critical energy density (today), in units of eV^4
		inline double DH () {return mDH;}          // Hubble distance in Mpc, c/H0
		inline double zCMB () {return mzcmb;}      // Redshift of recombination

	private:
		// Internal storage for values
		double mOmegaM;
		double mOmegaB;
		double mOmegaR;
		double mOmegaGamma;
		double mOmegaK;
		double mNeff;
		double mT;
		double mh;
		double mH0;
		double mz0;
		double mrhoc;
		double mDH;
		double mzcmb;

};

#endif /* PARAMS_H_ */
