/*
 * process.cpp
 *
 * Contains the routines to perform post-processing of data after the fact.
 *
 * At the moment, just computes the various distance measures used in cosmology.
 *
 */

#include "process.h"

using namespace std;

int PostProcessing(vector<double>& hubble, vector<double>& redshift, IntParams &params, Output &output) {
	// Takes in a vector of hubble values and z values (in ascending z order)
	// Also takes in some basic parameters about the cosmology
	// as well as the output class
	// Computes distance measures

	// Store the number of rows
	int numrows = hubble.size();

	// Step 1: Construct an interpolation of H and z
	// Trick: Get a pointer to the array for H and z
	double* H = &hubble[0];
	double* z = &redshift[0];

	// Initialize the spline
	splinetools myspline;
	myspline.acc = gsl_interp_accel_alloc();
	myspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);

	// Construct the spline
	gsl_spline_init (myspline.spline, z, H, numrows);

	// To get the value from the spline, use the following command
	// H(z) = gsl_spline_eval (myspline.spline, z, myspline.acc);


	// Step 2: Integrate
	// Now that we have our spline, it's time to integrate to obtain the appropriate quantities
	// We compute the quantity D_C(z)/D_H = \int_0^z dz'/H(z')/(1 + z')
	vector<double> DC;

	// Set up the integration stuff
	gsl_odeiv2_step *step;
	gsl_odeiv2_control *control;
	gsl_odeiv2_evolve *evolve;

	// Integration method as supplied by GSL
	const gsl_odeiv2_step_type *Type = gsl_odeiv2_step_rk8pd;

	// This is the number of elements that are being integrated
	size_t numelements = 1;

	// Initialize GSL
	step = gsl_odeiv2_step_alloc(Type, numelements);
	control = gsl_odeiv2_control_yp_new(0.0, 1e-12); // Only care about relative error (this is quite a stringent tolerance)
	evolve = gsl_odeiv2_evolve_alloc(numelements);

	// This is the stepsize that has been recommended by the integrator
	double stepsize;

	// Set the initial stepsize to be 10% of the first quantity
	stepsize = redshift[1] / 10;

	// Sets up the system to integrate, including the function that specifies the derivatives, NULL for the Jacobian, the number of elements
	// being integrated, and the parameters to pass through to the function
	gsl_odeiv2_system sys = { PPintfunc, NULL, numelements, &myspline };

	// We need to store the current value of the function we're integrating as an array
	double data[1] = {0.0};
	double currentz = 0.0;
	double zfinish;

	int status = GSL_SUCCESS;
	// We want to integrate from z = 0 to each z value in the data
	for (int i = 0; i < numrows; ++i) {
		// Read in the next z value to get to
		zfinish = redshift[i];

		// Integrate forwards in z until we get to the point we want
		while (currentz < zfinish) {

			// Do the integration
			status = gsl_odeiv2_evolve_apply(evolve, control, step, &sys, &currentz, zfinish, &stepsize, data);

			// Check for problems
			if (status != GSL_SUCCESS)
				break;

		}

		// Check for problems
		if (status != GSL_SUCCESS)
			break;

		// We're up to the next value of z; add it to the vector
		DC.push_back(data[0]);

	}

	// Check for problems
	if (status != GSL_SUCCESS) {
		std::stringstream erroroutput;
		erroroutput << "Error in integrating distance measures: " << status << std::endl;
		cout << erroroutput;
		output.printlog(erroroutput.str());
	}

	// Release the integrator memory
	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (step);

	// Release the spline memory
	gsl_spline_free (myspline.spline);
	gsl_interp_accel_free (myspline.acc);

	// Get out with an error response if necessary
	if (status != GSL_SUCCESS)
		return status;


	// Step 3: Compute the other distance measures
	// We have D_C(z)/D_H. The next step is to compute the rest of the measures (all / D_H, except for mu)
	vector<double> DM;
	vector<double> DA;
	vector<double> DL;
	vector<double> mu;
	// Reserve space in the vectors (because we can)
	DM.reserve(numrows);
	DA.reserve(numrows);
	DL.reserve(numrows);
	mu.reserve(numrows);

	// Extract OmegaK from the parameters
	double OmegaK = params.getparams().OmegaK();

	// The computation of DM depends on whether OmegaK is 0, positive or negative
	if (OmegaK > 0) {
		// k < 0
		double rootk = pow(OmegaK, 0.5);
		for (int i = 0; i < numrows; i++) {
			DM.push_back(sinh(rootk * DC[i]) / rootk);
		}
	}
	else if (OmegaK < 0) {
		// k > 0
		double rootk = pow((-OmegaK), 0.5);
		for (int i = 0; i < numrows; i++) {
			DM.push_back(sin(rootk * DC[i]) / rootk);
		}
	}
	else {
		// OmegaK == 0
		DM = DC;
	}

	// Next, compute DA and DL
	DA = DM;
	DL = DM;
	for (int i = 0; i < numrows; i++) {
		currentz = 1 + redshift[i];
		DA[i] /= currentz;
		DL[i] *= currentz;
	}

	// Finally, compute mu. This requires us to know h
	// In particular, we need to know c (in km/s) divided by h divided by 100 (ie, we want c/H0 in Mpc)
	double conh = 2997.92458 / params.getparams().h();
	for (int i = 0; i < numrows; i++) {
		mu.push_back(25 + 5 * log10 (DL[i] * conh));
	}


	// Step 4: Output all data. Convert all distance measurements to Mpc
	double DH = params.getparams().DH();
	for (int i = 0; i < numrows; i++) {
		output.postprintstep(redshift[i], hubble[i], DC[i] * DH, DM[i] * DH, DA[i] * DH, DL[i] * DH, mu[i]);
	}

	// Success!
	return 0;
}

int PPintfunc(double z, const double data[], double derivs[], void *params) {
	// This routine returns the derivative of the function we wish to integrate.
	// In this case, we're integrating \int_0^z dz'/H(z')/(1 + z')
	// Thus, the derivative is 1/H(z')/(1 + z')

	// Extract parameters
	splinetools myParams = *(splinetools *) params;

	// Calculate H(z)
	double H = gsl_spline_eval (myParams.spline, z, myParams.acc);

	// Calculate the derivative
	derivs[0] = 1 / (1.0 + z) / H;

	// Return success!
	return GSL_SUCCESS;

}
