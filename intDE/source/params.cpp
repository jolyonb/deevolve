/*
 * params.c
 *
 * Code supporting the reading and construction of cosmological parameters.
 *
 */

#include "params.h"

Parameters::Parameters(IniReader &init) {

	// A bunch of physical constants
	const double hubble0 = 3.24077929e-18; // This is 100 km/s/Mpc in units of s^-1
	const double GNewton = 6.67384e-11; // m^3 / kg / s^2
	const double c = 299792458; // m / s
	const double hbar = 1.05457173e-34; // Joule seconds
	const double pi = 3.14159265359;
	const double joulesinev = 6.24150934e18; // 1 joule in electron volts
	const double kinev = 8.61733238e-5; // Boltzmann's constant in eV/K

	// Read in a number of constants from the ini file
	mT = init.getiniDouble("Tgamma", 2.72548, "Cosmology");
	mz0 = init.getiniDouble("zInit", 1.0e4, "Cosmology");
	mNeff = init.getiniDouble("Neff", 3.046, "Cosmology");

    // Read in the cosmological values for matter and curvature
    mOmegaMh2 = init.getiniDouble("Omegamh2", 0.14305, "Cosmology");
    mOmegaBh2 = init.getiniDouble("Omegabh2", 0.022032, "Cosmology");
    mOmegaKh2 = init.getiniDouble("Omegakh2", 0.0, "Cosmology");
    // Make sure that fraction of baryonic matter is contained in fraction of matter
    if (mOmegaMh2 < mOmegaBh2) mOmegaBh2 = mOmegaMh2;

    // Compute the photon and radiation values
    // rho_0 = 3 H_0^2 / 8 pi G (in eV^4)
    double rho0 = 3.0 * pow(hubble0, 2.0) * pow(c, 5.0) * pow(hbar, 3.0) / 8 / pi / GNewton * pow(joulesinev, 4.0);
    // rho_\gamma = pi^2/15 * (kT)^4
    mOmegaGammah2 = pow(pi, 2.0) / 15 * pow(mT, 4.0) * pow(kinev, 4.0) / rho0;
	// But this is just for photons! For radiation, we need to add in neutrinos as well, which have the same energy density,
	// but at a lower temperature by a factor of (4/11)^(1/3), so total factor lower by
	// (4/11)^(4/3), but then increased by a factor of N_eff for the effective number of relativistic species
	// Also need a factor of 7/8 for fermions.
	mOmegaRh2 = mOmegaGammah2 * (1 + mNeff * 7.0 * pow(4.0/11.0, 4.0/3.0) / 8.0);

    // Store c / (100 km/s/Mpc) (in Mpc)
    mconh0 = c / 100000.0;

	// Set h to zero to indicate that it hasn't been initialized
	mh = 0.0;

}
