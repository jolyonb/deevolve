/*
 * main.cpp
 *
 * This is the program entry point for the program. The purpose of this software is to evolve scalar field models through
 * cosmological time.
 *
 * For library dependencies, see main.h
 *
 * Jolyon Bloomfield, March 2014
 *
 */

// All the includes are hidden in the header file
#include "main.h"

// Our program entry point
int main(int argc, char* argv[]) {

	// Read the params.ini file
	// If there is a command line argument, assume that it is a different file than params.ini
	IniReader inifile;
	if (argc > 1)
		inifile.read(argv[1]);
	else
		inifile.read("params.ini");

	// Set up the integrator
	Integrator *myIntegrator = new Integrator;

	// Set up the equations of motion/model class
//	Model *myModel = new Quintessence();
	Model *myModel = new QuintessenceH();
//	Model *myModel = new LambdaCDM();

	// Set up the parameters - OmegaM, Tgamma, OmegaK, z_init, h (of H_0 = h * 100 km/s/Mpc) and the model
	Parameters *myParams = new Parameters(inifile.getiniDouble("Omegam", 0.3, "general"),
			                              inifile.getiniDouble("Tgamma", 2.72548, "general"), 0.01, 1.0e4, 0.7);

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

	// The data array stores a, \phi, and \dot{phi} through the evolution. The fourth parameter is an initial value of \dot{a}/a,
	// which may be necessary in some models.
	double data[4] = { 1.0 / (1.0 + myParams->z0()), phi0, phidot0, 0.0 };

	// Set up the filenames to output
	string outputdir = "logs";
	string basename = "run";
	// Go and find our appropriate file name (using 4 digit numbers)
	string outputname = getfilename(outputdir, basename, 4);

	// Set up the output class
	Output *myOutput = new BasicDump(outputname);
	// Check that output is a go
	if (!myOutput->filesready()) {
		cerr << "Unable to open files for output." << endl;
		return -1;
	}

	// Set up the consistency check class
	// Consistency *myChecker = new Consistency();
	Consistency *myChecker = new SimpleCheck();

	// Allow the model to initialize itself
	int initresult = myModel->init(data, starttime, *myParams);
	if (initresult != 0) {
		cerr << "Unable to initialize model." << endl;
		return -1;
	}

	// Get the output class to write out information on the run
	myOutput->printinfo(data, *myIntParams);

	// Allow the model class to write out any information on the run
	myOutput->printlog(myModel->description());

	// Start timing!
	boost::timer::cpu_timer myTimer;

	// Do the evolution!
	int result = BeginEvolution(*myIntegrator, *myIntParams, data, starttime, endtime, *myOutput, *myChecker);

	// Print a goodbye message, using time in milliseconds
	myTimer.stop();
	myOutput->printfinish(myTimer.elapsed().wall / 1e6);

	// Clean up
	delete myChecker;
	delete myOutput;
	delete myIntParams;
	delete myParams;
	delete myModel;
	delete myIntegrator;

	// Exit gracefully with result from integration
	return result;
}

int BeginEvolution(Integrator &integrator, IntParams &params, double data[], double starttime, double endtime, Output &output, Consistency &check) {
	// This routine takes in a number of parameters, and performs the cosmological background evolution
	// integrator is the class that handles integration steps
	// params is the class that stores the cosmological parameters
	// data is an array of three elements that stores a, \phi, and \dot{\phi}. This is the data that gets evolved
	// starttime stores the starting value of conformal time for the evolution
	// endtime stores a value of conformal time after which we should stop the evolution
	// model is a class that computes the equations of motion

	// We need our own time variable to step forwards
	double time = starttime;
	// The result from the integrator. It returns GSL_SUCCESS (0) when everything works
	int result;
	// An array to hold the status information
	double status[17];

	// Write the column headers
	output.printheading(data, params);
	// Extract the initial state from the model
	params.getmodel().getstate(data, time, status, params.getparams());
	// Finally, write the initial conditions
	output.printstep(data, time, params, status);

	// This is the loop that takes successive integration steps. Halt if we go past the maximum evolution length.
	while (time < endtime) {

		// Take a step
		result = integrator.dointstep(intfunc, params, data, time, endtime);

		// If the step failed, break out of the loop
		if (result != GSL_SUCCESS)
			break;

		// Extract the state from the model
		params.getmodel().getstate(data, time, status, params.getparams());

		// Overwrite the hubble parameter from the status array
		// Either the status array was set to the hubble parameter already,
		// or the Hubble parameter is not being evolved, and was computed in status.
		data[3] = status[3];

		// Get the output class to write out the state of the system
		output.printstep(data, time, params, status);

		// Take a look at the consistency of the data
		check.checkstate(data, time, params, output, status);

		// If we've shot past a = 1, then get out
		if (data[0] > 1.0)
			break;
	}

	// Take a look at the consistency of the final state
	check.checkfinal(data, time, params, output, status);

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

string getfilename(std::string &dir, std::string &filebase, int padding) {
	// This routine takes in a directory and an output name
	// It goes and finds the first available filename of the form dir / filebase 00001 etc
	// eg., dir/run00001.log and dir/run00001.dat
	// It checks that both files are free
	// Note that even if the output is going to screen, this routine won't make anything bad happen

	using namespace boost::filesystem;

	// Firstly, make sure that the directory exists
	if (!exists(dir + "/")) {
		// Directory doesn't exist. Make it.
		create_directory(dir);
		std::cout << "Creating directory " << dir << "/" << endl;
	}

	// Secondly, find a unique filename
	for (int counter = 1; ; counter++) {
		// Construct the file number
		string filenum;
		ostringstream convert;
		convert << counter;
		filenum = convert.str();
		// Pad the file number with the appropriate number of zeroes
		int len = filenum.length();
		for (int i = 0; i < padding - len; i++)
			filenum = "0" + filenum;

		// Check for the files
		if (exists(dir + "/" + filebase + filenum + ".log"))
			continue;
		if (exists(dir + "/" + filebase + filenum + ".dat"))
			continue;

		// If we got to here, we have a unique filename; return it
		return dir + "/" + filebase + filenum;
	}

	// We really shouldn't get here, but parsers like making sure there's a return
	return "error";

}
