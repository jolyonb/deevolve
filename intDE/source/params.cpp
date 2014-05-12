/*
 * params.c
 *
 * Code supporting the reading and construction of cosmological parameters.
 *
 */

#include "params.h"

Parameters::Parameters(IniReader &init) {

	// A bunch of physical constants
	const double hubble0conv = 3.24077929e-18; // This is 100 km/s/Mpc in units of s^-1
	const double GNewton = 6.67384e-11; // m^3 / kg / s^2
	const double c = 299792458; // m / s
	const double hbar = 1.05457173e-34; // Joule seconds
	const double pi = 3.14159265359;
	const double joulesinev = 6.24150934e18; // 1 joule in electron volts
	const double kinev = 8.61733238e-5; // Boltzmann's constant in eV/K

	// Read in a number of constants from the ini file
	mOmegaK = init.getiniDouble("Omegak", 0.0, "Cosmology");
	mT = init.getiniDouble("Tgamma", 2.72548, "Cosmology");
	mh = init.getiniDouble("Hubbleh", 0.7, "Cosmology");
	double mh2 = mh * mh;
	mz0 = init.getiniDouble("zInit", 1.0e4, "Cosmology");
	mH0 = hubble0conv * mh;
	mNeff = init.getiniDouble("Neff", 3.046, "Cosmology");

	// Do we read in OmegaB * h^2 or just OmegaB? (and same for OmegaM)
	double mOmegaMh2;
	double mOmegaBh2;
	if (init.getiniBool("Useh2", false, "Cosmology")) {
		mOmegaMh2 = init.getiniDouble("Omegamh2", 0.14305, "Cosmology");
		mOmegaBh2 = init.getiniDouble("Omegabh2", 0.022032, "Cosmology");
		mOmegaM = mOmegaMh2 / mh2;
		mOmegaB = mOmegaBh2 / mh2;
	} else {
		mOmegaM = init.getiniDouble("Omegam", 0.3, "Cosmology");
		mOmegaB = init.getiniDouble("Omegab", 0.05, "Cosmology");
		mOmegaMh2 = mOmegaM * mh2;
		mOmegaBh2 = mOmegaB * mh2;
	}
	if (mOmegaM < mOmegaB) {
		// Make sure that fraction of baryonic matter is contained in fraction of matter
		mOmegaB = mOmegaM;
		mOmegaBh2 = mOmegaMh2;
	}

	// The quantity 3 H_0^2 / 8 pi G = 3 H_0^2 c^5 hbar^3 / 8 pi G has units of Joules^4 in the above units
	// This can be converted into eV^4, which is the units we will express the critical density in
	mrhoc = 3.0 * pow(mH0, 2.0) * pow(c, 5.0) * pow(hbar, 3.0) / 8 / pi / GNewton * pow(joulesinev, 4.0);

	// rho_\gamma = pi^2/15 * T^4
	mOmegaGamma = pow(pi, 2.0) / 15 * pow(mT, 4.0) / mrhoc * pow(kinev, 4.0);
	// But this is just for photons! For radiation, we need to add in neutrinos as well, which have the same energy density,
	// but at a lower temperature by a factor of (4/11)^(1/3), so total factor lower by
	// (4/11)^(4/3), but then increased by a factor of 3 (three species of neutrino;
	// v and vbar, but photon already has two polarizations included)
	// Also need a factor of 7/8 for fermions.
	mOmegaR = mOmegaGamma * (1 + mNeff * 7.0 * pow(4.0/11.0, 4.0/3.0) / 8.0);

	// Calculate DH = c/H0 in Mpc (used in distance measurements)
	mDH = c / 100000.0 / mh;

	// Calculate redshift of recombination (CMB formation)
	// Computed using Eq E.1 from astro-ph/9510117v2 (note that Omega_0 = Omega_m) for percent level accuracy
	double g1;
	double g2;
	double OmegaBh2 = mOmegaB * mh2;
	g1 = 0.0783 * pow( OmegaBh2, -0.238) / (1.0 + 39.5 * pow(OmegaBh2, 0.763));
	g2 = 0.560 / (1.0 + 21.1 * pow(OmegaBh2, 1.81));
	mzcmb = 1048 * (1.0 + 0.00124 * pow(OmegaBh2, -0.738)) * (1.0 + g1 * pow(mOmegaM * mh2, g2));

}
