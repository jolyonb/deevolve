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
#include <vector>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

// We need a data type to pass a spline as well as the spline accelerator into the integration function
struct splinetools {
	gsl_interp_accel *acc;
	gsl_spline *spline;
};

// Routine that performs all postprocessing
int PostProcessing(vector<double>&, vector<double>&, IntParams&, Output&);

// Routine that returns the integrand
int PPintfunc(double, const double*, double*, void*);

#endif /* PROCESS_H_ */
