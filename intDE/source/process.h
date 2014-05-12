/*
 * process.h
 *
 * Describes the routines used for performing post-processing
 */

#ifndef PROCESS_H_
#define PROCESS_H_

#include "intparams.h"
#include "output.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <boost/filesystem.hpp> // Used for making sure sn1a file exists
#include <iomanip>


// We need a data type to pass a spline as well as the spline accelerator into the integration function
struct splinetools {
	gsl_interp_accel *acc;
	gsl_spline *spline;
	double param; // Also allow an extra parameter to go through to the integrator
};

// Routine that performs all postprocessing
int PostProcessingDist(vector<double>&, vector<double>&, vector<double>&, vector<double>&,
		vector<double>&, vector<double>&, vector<double>&, double &rs, IntParams&, Output&, IniReader &init);

// Routines that return the integrand for various integrals
int PPintfunc(double, const double*, double*, void*);
double rsintfunc(double z, void *params);
double rsintfuncinf(double z, void *params);

// Routine to compute chi^2 values for supernovae data
int chi2SN1a(vector<double>& redshift, vector<double>& mu, Output &output, IniReader &init);

// Routine to compute chi^2 values for CMB data
int chi2CMB(vector<double>& redshift, vector<double>& mu, double &rs, Output &output, IntParams &params);

// Routine to compute chi^2 values for BAO data
//int chi2BAO(vector<double>& redshift, vector<double>& mu, Output &output);

// Routine to compute chi squared values
int calcchisquared(vector<double>& redshift,
		           vector<double>& hubble,
		           vector<double>& DC,
		           vector<double>& DM,
		           vector<double>& DA,
		           vector<double>& DL,
		           vector<double>& mu,
		           double rs,
		           IntParams &params, Output &output, IniReader &init);

// Routines to compute chi^2 for WMAP and Planck distance posteriors
double chi2WMAP (double lA, double R, double z);
double chi2Planck (double lA, double R, double z);

#endif /* PROCESS_H_ */
