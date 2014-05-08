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
#include <boost/filesystem.hpp> // Used for making sure sn1a file exists


// We need a data type to pass a spline as well as the spline accelerator into the integration function
struct splinetools {
	gsl_interp_accel *acc;
	gsl_spline *spline;
};

// Routine that performs all postprocessing
int PostProcessing(vector<double>&, vector<double>&, IntParams&, Output&, IniReader &init);

// Routine that returns the integrand
int PPintfunc(double, const double*, double*, void*);

// Routine to compute chi squared values
int calcchisquared(vector<double>& redshift,
		           vector<double>& hubble,
		           vector<double>& DC,
		           vector<double>& DM,
		           vector<double>& DA,
		           vector<double>& DL,
		           vector<double>& mu,
		           IntParams &params, Output &output, IniReader &init);

#endif /* PROCESS_H_ */
