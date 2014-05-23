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
using namespace boost::filesystem;

int PostProcessingDist(vector<double>& hubble,
		vector<double>& redshift,
		vector<double>& DC,
		vector<double>& DM,
		vector<double>& DA,
		vector<double>& DL,
		vector<double>& mu,
		double &rs,
		IntParams &params, Output &output, IniReader &init) {
	// Takes in a vector of hubble values and z values (in ascending z order)
	// The rest of the vectors are assumed to be empty
	// rs will be given the sound horizon scale at z_CMB
	// Also takes in some basic parameters about the cosmology
	// as well as the output class
	// Computes distance measures
	// If requested in the ini file, also computes chi squared values

	// Store the number of rows
	int numrows = hubble.size();

	// Step 1: Construct an interpolation of H and z
	// Trick: Get a pointer to the array for H and z
	double* pH = &hubble[0];
	double* pz = &redshift[0];

	// Initialize the spline
	splinetools myspline;
	myspline.acc = gsl_interp_accel_alloc();
	myspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);

	// Construct the spline
	gsl_spline_init (myspline.spline, pz, pH, numrows);

	// To get the value from the spline, use the following command
	// H(z) = gsl_spline_eval (myspline.spline, z, myspline.acc);


	// Step 2: Integrate
	// Now that we have our spline, it's time to integrate to obtain the appropriate quantities
	// We compute the quantity D_C(z)/D_H = \int_0^z dz'/H(z')/(1 + z')

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
	stepsize = redshift[1] / 10.0;

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

	// Release the integrator memory
	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (step);

	// Get out with an error response if necessary
	if (status != GSL_SUCCESS) {
		// Error message
		std::stringstream erroroutput;
		erroroutput << "Error in integrating distance measures: " << status << std::endl;
		cout << erroroutput;
		output.printlog(erroroutput.str());

		// Release the spline memory
		gsl_spline_free (myspline.spline);
		gsl_interp_accel_free (myspline.acc);
		return status;
    }


	// Step 3: Compute the other distance measures
	// We have D_C(z)/D_H. The next step is to compute the rest of the measures (all / D_H, except for mu)

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

	// At z = 0, all distance measures are zero, H = 1, and mu = -infinity.
	// This is a problem for mu, which becomes infinite.
	// For this reason, we remove the first entry of each distance measurement when reporting data.

	// Step 4: Output all data. Convert all distance measurements to Mpc
	double DH = params.getparams().DH();
	for (int i = 1; i < numrows; i++) {
		output.postprintstep(redshift[i], hubble[i], DC[i] * DH, DM[i] * DH, DA[i] * DH, DL[i] * DH, mu[i]);
	}


	// Step 5: Calculate the sound horizon distance at recombination
	// This one is a little more complicated.
	// The expression we need is the following.
	// r_s / D_H = \frac{1}{\sqrt{3}} \int^{\infty}_{z_{CMB}} \frac{dz'}{(1 + z') \tilde{\ch}(z')} \frac{1}{\sqrt{1 + 3 \Omega_B / 4 \Omega_\gamma (1 + z')}}
	// Start by computing the quantity that depends on OmegaB and OmegaR, and storing it in the integration helper
	double alpha = 3 * params.getparams().OmegaB() / 4 / params.getparams().OmegaGamma();
	myspline.param = alpha;
	// Note that this integral goes from zCMB to infinity. Our evolution doesn't go to infinite redshift however.
	// We need to evaluate the integral over our span of redshift, as well as over the infinite portion.
	// As we don't have information outside our span, we use a LambdaCDM prediction for that period, which should be accurate for most purposes.

	// As we are only computing one value here, rather than a value at various redshifts, we can use a specialized integration routine.
	// Set up the integration workspace
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000);

	// First, we integrate over the data that we have
	double intresult, error;
	gsl_function intFunc;
	intFunc.function = &rsintfunc;
	intFunc.params = &myspline;

	// Perform the integration
	gsl_integration_qag (&intFunc, params.getparams().zCMB(), params.getparams().z0(), 1e-8, 0, 1000, 6, workspace, &intresult, &error);
	// Store the result
	rs = intresult;

	// Next, we integrate over the infinite part of the integral. Note that we need to pass in all the parameters to calculate this integrand.
	intFunc.function = &rsintfuncinf;
	intFunc.params = &params;

	// Perform the integration in two steps. Once, to redshift 1e5, using the more accurate integrator. Then, to infinite redshift.
	// First, finite part
	gsl_integration_qag (&intFunc, params.getparams().z0(), 1e5, 1e-8, 0, 1000, 6, workspace, &intresult, &error);
	// Store the result
	rs += intresult;
	// Infinite part.
	gsl_integration_qagiu (&intFunc, 1e5, 1e-8, 0, 1000, workspace, &intresult, &error);
	// Store the result
	rs += intresult;

	// Include the multiplicative factor of 1/sqrt(3)
	rs *= pow(3.0, -0.5);

	// Release the memory from the integrator
	gsl_integration_workspace_free (workspace);

	// Release the spline memory
	gsl_spline_free (myspline.spline);
	gsl_interp_accel_free (myspline.acc);

	// Report the result
	{
		std::stringstream reportoutput;
		reportoutput << setprecision(8) << "Sound horizon scale at recombination: " << rs * DH << " Mpc" << std::endl;
		output.printlog(reportoutput.str());
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

double rsintfunc(double z, void *params) {
	// This routine returns the derivative of the function we wish to integrate.
	// In this case, we're integrating \int dz / (1 + z') H(z') sqrt{1 + alpha / (1 + z')}
	// Thus, the derivative is 1 / (1 + z') H(z') sqrt{1 + alpha / (1 + z')}

	// Extract parameters
	splinetools myParams = *(splinetools *) params;

	// Calculate H(z)
	double H = gsl_spline_eval (myParams.spline, z, myParams.acc);

	// Calculate the derivative
	return 1 / (1.0 + z) / H / pow(1 + myParams.param / (1.0 + z), 0.5);

}

double rsintfuncinf(double z, void *params) {
	// This routine returns the derivative of the function we wish to integrate.
	// In this case, we're integrating \int dz / (1 + z') H(z') sqrt{1 + alpha / (1 + z')}
	// Thus, the derivative is 1 / (1 + z') H(z') sqrt{1 + alpha / (1 + z')}
	// However, here we do not have the benefit of a spline

	// Extract parameters
	IntParams myParams = *(IntParams *) params;

	// Calculate H(z) based on LambdaCDM FRW evolution with radiation, matter and curvature
	double a = 1.0 + z;
	double H = pow(a * a * myParams.getparams().OmegaR() + a * myParams.getparams().OmegaM() + myParams.getparams().OmegaK(), 0.5);
	double alpha = 3 * myParams.getparams().OmegaB() / 4 / myParams.getparams().OmegaR();

	// Calculate the derivative
	return 1 / a / H / pow(1 + alpha / a, 0.5);

}

int chi2SN1a(vector<double>& redshift, vector<double>& mu, Output &output, IniReader &init) {
	// This routine computes the chi^2 value for supernovae measurements

	// Step 1: We will need to know mu at specific redshifts. Construct the required interpolater.

	// Trick: Get pointers to the various arrays
	// At z = 0, mu = -infinity.
	// For this reason, we remove the first entry when constructing the interpolation.
	double* pz = &redshift[1];
	double* pmu = &mu[1];
	int numrows = redshift.size() - 1;

	// Initialize the splines
	// Distance modulus is used for the SN1a data
	splinetools muspline;
	muspline.acc = gsl_interp_accel_alloc();
	muspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);

	// Construct the spline
	gsl_spline_init (muspline.spline, pz, pmu, numrows);

	// Step 2: Read in all of the data from the SN1a data
	string sn1afile = init.getiniString("union21", "SCPUnion2.1_mu_vs_z.txt", "Function");
	// Check that the file exists before further processing
	if (exists(sn1afile)) {
		// Begin by reading in the data file
		ifstream f(sn1afile.c_str());
		string l;
		vector<vector<double> > rows; // This is a vector of vectors

		// Read in the file one line at a time into the string l
		while(getline(f, l)) {
			// If the first character is a # (comment), move onto the next line
			if (l[0] == '#')
				continue;
			// Convert the string l into a stringstream s
			stringstream s(l);
			string extract;
			double entry;
			vector<double> row;
			// Extract entries one at a time, using a tab as a delimiter
			// Ignore the first entry, which is the supernovae name
			bool first = true;
			while(getline(s, extract, '\t')) {
				if (first) {
					first = false;
				} else {
				    row.push_back(atof(extract.c_str()));
				}
			}
			rows.push_back(row);
		}

		// rows now contains a vector of doubles. Each row is the data from a supernovae.
		// Each row contains four pieces of information: redshift, distance modulus, error on distance modulus, and a piece of
		// information that is not of use to us.

		// Iterate through each row, and construct the chi^2
		double chi2 = 0;
		double val = 0;
		int numrows = rows.size();
		for (int i = 0; numrows; i++) {
			// Add (mu - mu(z))^2 / sigma^2 to the chi^2
			val = (rows[i][1] - gsl_spline_eval (muspline.spline, rows[i][0], muspline.acc)) / rows[i][2];
			chi2 += val * val;
		}

		// Having gotten here, present the results
		std::stringstream printing;
		printing << "SN1a chi^2 = " << chi2 << " (computed from " << rows.size() << " data points)" << endl;
		output.printlog(printing.str());
		// The minimum chi^2 for LambdaCDM for this data set is ~562.5

	}
	else {
		// Could not find data file
		std::stringstream printing;
		printing << "Error: cannot find Union2.1 SN1a data file." << endl;
		cout << printing.str();
		output.printlog(printing.str());
	}

	// Release the spline memory
	gsl_spline_free (muspline.spline);
	gsl_interp_accel_free (muspline.acc);

	return 0;
}

int chi2CMB(vector<double>& redshift, vector<double>& DA, double &rs, Output &output, IntParams &params) {
	// This routine computes the chi^2 value for CMB distance posteriors

	// We use the WMAP distance posteriors on the acoustic scale and the shift parameter
	// See http://lambda.gsfc.nasa.gov/product/map/dr3/pub_papers/fiveyear/cosmology/wmap_5yr_cosmo_reprint.pdf
	// This needs the angular diameter distance at recombination, so we'll need a new interpolater
	// Extract the appropriate data
	double* pz = &redshift[0];
	double* pDA = &DA[0];
	int numrows = redshift.size();

	// Construct the spline
	splinetools DAspline;
	DAspline.acc = gsl_interp_accel_alloc();
	DAspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);
	gsl_spline_init (DAspline.spline, pz, pDA, numrows);

	// Get the redshift of recombination
	double zCMB = params.getparams().zCMB();

	// Calculate l_A (note that as we're taking the ratio between the two distance scales, D_H doesn't appear here at all)
	double la = (1.0 + zCMB) * gsl_spline_eval (DAspline.spline, zCMB, DAspline.acc) * M_PI / rs;

	// Calculate R (note that we don't need to divide by D_H here)
	double R = pow(params.getparams().OmegaM(), 0.5) * (1.0 + zCMB) * gsl_spline_eval (DAspline.spline, zCMB, DAspline.acc);

	// Release the spline memory
	gsl_spline_free (DAspline.spline);
	gsl_interp_accel_free (DAspline.acc);

	// Finally, compute the chi^2 values for WMAP and Planck
	double chi2wmap = chi2WMAP(la, R, zCMB);
	double chi2p = chi2Planck(la, R, zCMB);

	// Output everything
	{
		std::stringstream printing;
		printing << setprecision(8) << "Acoustic scale: l_A = " << la << endl;
		printing << "Shift parameter: R = " << R << endl;
		printing << "Redshift of CMB: z_CMB = " << zCMB << endl << endl;

		printing << "WMAP Distance Posterior chi^2 = " << setprecision(8) << chi2wmap << endl;
		printing << "Planck Distance Posterior chi^2 = " << setprecision(8) << chi2p << endl;
		output.printlog(printing.str());
	}

	// Return success
	return 0;

}

// Computes the chi^2 value for the CMB distance posteriors from WMAP 9 data
// See Table 11 in arxiv.org/pdf/1212.5226v3.pdf and numbers below
double chi2WMAP (double lA, double R, double z) {
	// Construct deltas
	double deltal = lA - 302.4;
	double deltaR = R - 1.7246;
	double deltaz = z - 1090.88;

	double chi2 = 0;
	// Add contributions from diagonals first
	chi2 += deltal * deltal * 3.182;
	chi2 += deltaR * deltaR * 11887.879;
	chi2 += deltaz * deltaz * 4.556;
	// Add contributions from off-diagonal terms next
	chi2 += 2 * deltal * deltaR * 18.253;
	chi2 += - 2 * deltal * deltaz * 1.429;
	chi2 += - 2 * deltaR * deltaz * 193.808;

	// Return the result
	return chi2;
}

// Computes the chi^2 value for the CMB distance posteriors from Planck 1 data
// See Tables 1 and 2 in http://arxiv.org/pdf/1309.0679v1.pdf
// Note, we use the LambdaCDM values, as specified for the covariance matrix in the text
double chi2Planck (double lA, double R, double z) {
	// Construct deltas
	double deltal = lA - 301.77;
	double deltaR = R - 1.7477;
	double deltaz = z - 1090.25;

	double chi2 = 0;
	// Add contributions from diagonals first
	chi2 += deltal * deltal * 44.077;
	chi2 += deltaR * deltaR * 48976.330;
	chi2 += deltaz * deltaz * 12.592;
	// Add contributions from off-diagonal terms next
	chi2 += - 2 * deltal * deltaR * 383.927;
	chi2 += - 2 * deltal * deltaz * 1.941;
	chi2 += - 2 * deltaR * deltaz * 630.791;

	// Return the result
	return chi2;
}
