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
#include <iostream>

class Parameters {
	public:

		// Constructor. Sets the parameters of the model, which cannot be changed afterwards
		Parameters(const double OmegaM, const double OmegaB, const double Tgamma, const double OmegaK, const double z0, const double h) {

			const double hubble0conv = 3.24e-18; // This is 100 km/s/Mpc in units of s^-1

			mOmegaM = OmegaM;
			mOmegaB = OmegaB;
			if (mOmegaM < mOmegaB) // Make sure that fraction of baryonic matter is contained in fraction of matter
				mOmegaB = mOmegaM;
			mOmegaK = OmegaK;
			mT = Tgamma;
			mh = h;
			mz0 = z0;
			mH0 = hubble0conv * h;

			const double GNewton = 6.67384e-11; // m^3 / kg / s^2
			const double c = 299792458; // m / s
			const double hbar = 1.05457173e-34; // Joule seconds
			const double pi = 3.14159265359;
			const double joulesinev = 6.24150934e18; // 1 joule in electron volts
			// The quantity 3 H_0^2 / 8 pi G = 3 H_0^2 c^5 hbar^3 / 8 pi G has units of Joules^4 in the above units
			// This can be converted into eV^4, which is the units we will express the critical density in
			mrhoc = 3.0 * pow(mH0, 2.0) * pow(c, 5.0) * pow(hbar, 3.0) / 8 / pi / GNewton * pow(joulesinev, 4.0);

			// rho_\gamma = pi^2/15 * T^4
			const double kinev = 8.61733238e-5; // eV/K
			mOmegaR = pow(pi, 2.0) / 15 * pow(Tgamma, 4.0) / mrhoc * pow(kinev, 4.0);
			// But this is just for photons! We need to add in neutrinos as well, which have the same energy density,
			// but at a lower temperature by a factor of (4/11)^(1/3), so total factor lower by
			// (4/11)^(4/3), but then increased by a factor of 3 (three species of neutrino;
			// v and vbar, but photon already has two polarizations included)
			// Also need a factor of 7/8 for fermions.
			mOmegaR *= 1 + 21.0 * pow(4.0/11.0, 4.0/3.0) / 8.0;

			// Calculate DH = c/H0 in Mpc (used in distance measurements)
			mDH = c / 1000.0 / 100.0 / mh;

			// Calculate redshift of recombination (CMB formation)
			// Computed using Eq E.1 from astro-ph/9510117v2 (note that Omega_0 = Omega_m) for percent level accuracy
			double g1;
			double g2;
			double OmegaBh2 = mOmegaB * mh * mh;
			g1 = 0.0783 * pow( OmegaBh2, -0.238) / (1 + 39.5 * pow(OmegaBh2, 0.763));
			g2 = 0.560 / (1 + 21.1 * pow(OmegaBh2, 1.81));
			mzcmb = 1048 * (1 + 0.00124 * pow(OmegaBh2, -0.738)) * (1 + g1 * pow(mOmegaM * mh * mh, g2));

		}

		// Getters for the parameters of the model
		inline double OmegaM () {return mOmegaM;}  // Energy fraction of matter today
		inline double OmegaB () {return mOmegaB;}  // Energy fraction of baryonic matter today
		inline double OmegaR () {return mOmegaR;}  // Energy fraction of radiation today
		inline double OmegaK () {return mOmegaK;}  // Energy fraction of spatial curvature today
		inline double Tgamma () {return mT;}       // Temperature of photons today
		inline double h () {return mh;}            // Hubble parameter is h * 100 km/s/Mpc
		inline double H0 () {return mH0;}          // Hubble parameter today (in s^-1)
		inline double z0 () {return mz0;}          // Redshift to start the evolution at
		inline double rhoc () {return mrhoc;}      // Critical energy density (today), in units of eV^4
		inline double DH () {return mDH;}          // Hubble distance in Mpc, c/H0
		inline double zCMB () {return mzcmb;}      // Redshift of recombination

	private:
		// Internal storage for values
		double mOmegaM;
		double mOmegaB;
		double mOmegaR;
		double mOmegaK;
		double mT;
		double mh;
		double mH0;
		double mz0;
		double mrhoc;
		double mDH;
		double mzcmb;

};

#endif /* PARAMS_H_ */
