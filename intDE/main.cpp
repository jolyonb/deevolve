/*
 * main.cpp
 *
 * This is the program entry point for the program. The purpose of this software is to evolve scalar field models through
 * cosmological time.
 *
 * Jolyon Bloomfield, March 2014
 *
 */

// All the includes are hidden in the header file
#include "main.h"

// Our program entry point. Ignore command line arguments.
int main(void) {

	// Set up the integrator
	Integrator *myIntegrator = new Integrator;

	// Set up the equations of motion/model class
	Model *myModel = new Quintessence();
//	Model *myModel = new LambdaCDM();

	// Set up the parameters - OmegaM, Tgamma, OmegaK, z_init, h (of H_0 = h * 100 km/s/Mpc) and the model
	Parameters *myParams = new Parameters(0.3, 2.72548, 0.01, 1.0e4, 0.7);

	// Load the model and parameters into a class to pass into the integration routine
	IntParams *myIntParams = new IntParams(*myParams, *myModel);

	// Set up the initial conditions
	// We use conformal time for our time coordinate. To make things easier, we start at tau = 0.
	double starttime = 0.0;
	// We don't want to integrate forever; specify an endtime in conformal time. I've chosen 100 here; this should
	// be chosen to correspond to about twice as long as you'd except the integration to take.
	// Note that we'll actually be checking to stop integration when the energy fraction of matter is roughly correct (a >= 1)
	double endtime = 10.0;
	// The initial conditions is specified by values for the initial a, \phi, and \dot{\phi}
	// where the overdot is a derivative with respect to conformal time.
	// Rather than setting a, it's simpler to set a redshift to begin at, and compute a from there, noting that 1 + z = a_{today}/a
	// and a_{today} = 1 by choice. This choice is specified in myParams.
	// \phi and \dot{\phi} will need to be set based on the model
	double phi0 = 0.01;
	double phidot0 = 0.0;

	// The data array stores a, \phi, and \dot{phi} through the evolution
	double data[3] = { 1.0 / (1.0 + myParams->z0()), phi0, phidot0 };

	// Set up the output class
	Output *myOutput = new BasicDump();

	// Do the evolution!
	int result = BeginEvolution(*myIntegrator, *myIntParams, data, starttime, endtime, *myOutput);

	// Clean up
	delete myOutput;
	delete myIntParams;
	delete myParams;
	delete myModel;
	delete myIntegrator;

	// Exit gracefully with result from integration
	return result;
}

int BeginEvolution(Integrator &integrator, IntParams &params, double data[], double starttime, double endtime, Output &output) {
	// This routine takes in a number of parameters, and performs the cosmological background evolution
	// integrator is the class that handles integration steps
	// params is the class that stores the cosmological parameters
	// data is an array of three elements that stores a, \phi, and \dot{\phi}. This is the data that gets evolved
	// starttime stores the starting value of conformal time for the evolution
	// endtime stores a value of conformal time after which we should stop the evolution
	// model is a class that computes the equations of motion

	// We need our own time variable to step forwards
	double time = starttime;
	// The result from the integrator. It returns GSL_SUCCESS (aka 0) when everything works
	int result;

	// Get the output class to write out information on the run
	output.printinfo(data, params);
	// Then, write the column headers
	output.printheading(data, params);
	// Finally, write the initial conditions
	output.printstep(data, time, params);

	// This is the loop that takes successive integration steps. Halt if we go past the maximum evolution length.
	while (time < endtime) {

		// Take a step
		result = integrator.dointstep(intfunc, params, data, time, endtime);

		// If the step failed, break out of the loop
		if (result != GSL_SUCCESS)
			break;

		// Get the output class to write out the state of the system
		output.printstep(data, time, params);

		// If we've shot past a = 1, then get out
		if (data[0] > 1.0)
			break;
	}

	return result;
}

int intfunc(double t, const double data[], double derivs[], void *params) {
	// This routine calculates the derivatives for a, phi and \dot{\phi}
	// It is called by the integration routine.

	// Extract parameters
	IntParams myParams = *(IntParams *) params;

	// Call the derivatives routine in the model to calculate the derivatives appropriately
	// Note that they don't depend on time
	int result = myParams.getmodel().derivatives(data, derivs, myParams.getparams());

	// Return the result
	return result;
}
