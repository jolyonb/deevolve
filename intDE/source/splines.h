/*
 * Stores information on how to make splines. Used in a few places.
 */

#ifndef SPLINES_H_
#define SPLINES_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// We need a data type to pass a spline as well as the spline accelerator into the integration function
struct splinetools {
	gsl_interp_accel *acc;
	gsl_spline *spline;
	double param; // Also allow an extra parameter to go through to the integrator
};

#endif /* SPLINES_H_ */
