/*
 * main.cpp
 *
 * This is the program entry point for the program. The purpose of this software is to evolve scalar field models through
 * cosmological time.
 *
 * For library dependencies, see main.h
 *
 * For info on rameters, see params.ini
 *
 * Most initialization (and memory release) is handled in the "main" routine.
 * All integration is handled in the "BeginEvolution" routine.
 * intfunc is a helper function for integration.
 * getfilename is a helper function for creating log files.
 *
 * Jolyon Bloomfield, March 2014
 *
 */

// All the includes are hidden in the header file
#include "main.h"

// Our program entry point
int main(int argc, char* argv[]) {

	// Just a number for return values
	int result;

	// Read the params.ini file
	// If there is a command line argument, assume that it is a different file than params.ini
	IniReader inifile;
	if (argc > 1)
		inifile.read(argv[1]);
	else
		inifile.read("params.ini");

	// Set up the integrator
	Integrator *myIntegrator = new Integrator;

	// Set up the model class
	Model *myModel;
	std::string parsestring = inifile.getiniString("model", "LambdaCDM", "Cosmology");
	if (parsestring == "Quintessence")
		myModel = new Quintessence();
	else
		myModel = new LambdaCDM();    // LambdaCDM is the default

	// Set up the parameters - OmegaM, Tgamma, OmegaK, z_init and h (of H_0 = h * 100 km/s/Mpc)
	Parameters *myParams = new Parameters(inifile.getiniDouble("Omegam", 0.3, "Cosmology"),
			                              inifile.getiniDouble("Tgamma", 2.72548, "Cosmology"),
			                              inifile.getiniDouble("Omegak", 0.01, "Cosmology"),
			                              inifile.getiniDouble("zInit", 1.0e4, "Cosmology"),
			                              inifile.getiniDouble("Hubbleh", 0.7, "Cosmology"));

	// Load the model and parameters into a class to pass into the integration routine
	IntParams *myIntParams = new IntParams(*myParams, *myModel);

	// Set up more initial conditions
	double starttime = inifile.getiniDouble("starttime", 0.0, "Function");
	double endtime = starttime + inifile.getiniDouble("maxtime", 10.0, "Function");
	double phi0 = inifile.getiniDouble("phi0", 0.01, "Cosmology");
	double phidot0 = inifile.getiniDouble("phidot0", 0.0, "Cosmology");;

	// The data array stores a, \phi, \dot{phi} and H through the evolution
	// H is calculated in the initialization of the model
	// The initial value of a is extracted from the starting redshift
	double data[4] = { 1.0 / (1.0 + myParams->z0()), phi0, phidot0, 0.0 };

	// Set up the filenames to output
	string outputdir = inifile.getiniString("logdir", "logs", "Function");
	string basename = inifile.getiniString("runname", "run", "Function");
	// Go and find our appropriate file name (using 4 digit numbers)
	string outputname = getfilename(outputdir, basename, inifile.getiniInt("numberpad", 4, "Function"));

	// Set up the output class
	Output *myOutput = new BasicDump(outputname);
	// Check that output is a go
	if (!myOutput->filesready()) {
		cerr << "Unable to open files for output." << endl;
		return -1;
	}

	// Set up the consistency check class
	parsestring = inifile.getiniString("consistencyclass", "None", "Function");
	Consistency *myChecker;
	if (parsestring == "SimpleCheck")
		myChecker = new SimpleCheck();
	else
		myChecker = new Consistency();  // Default option, which has no checking

	// Get the output class to write out information on the run
	myOutput->printinfo(data, *myIntParams);

	// Allow the model to initialize itself
	std::string response = myModel->init(data, starttime, *myParams, inifile, result);
	myOutput->printlog(response);
	if (result != 0) {
		cerr << "Unable to initialize model." << endl << response << endl;
		return -1;
	}

	// Start timing!
	boost::timer::cpu_timer myTimer;

	// Do the evolution!
	result = BeginEvolution(*myIntegrator, *myIntParams, data, starttime, endtime, *myOutput, *myChecker);

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

int BeginEvolution(Integrator &integrator, IntParams &params, double data[], const double starttime, const double endtime, Output &output, Consistency &check) {
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
	// And a double to hold the stepsize
	double stepsize;

	// Write the column headers
	output.printheading();
	// Extract the initial state from the model
	params.getmodel().getstate(data, time, status, params.getparams());
	// Write the initial conditions
	output.printstep(data, time, params, status);
	// Do a consistency check on the initial conditions
	check.checkstate(data, time, params, output, status);

	// This is the loop that takes successive integration steps. Halt if we go past the maximum evolution length.
	while (time < endtime) {

		// Take a step
		result = integrator.dointstep(intfunc, params, data, time, endtime);

		// If the step failed, break out of the loop
		if (result != GSL_SUCCESS) {
			cerr << "Integration routine failed." << endl;
			output.printlog("Integration routine failed.");
			break;
		}

		// Extract the state from the model
		params.getmodel().getstate(data, time, status, params.getparams());

		// Get the output class to write out the state of the system
		output.printstep(data, time, params, status);

		// Take a look at the consistency of the data
		check.checkstate(data, time, params, output, status);

		// If we've shot past a = 1, then get out
		if (data[0] > 1.0)
			break;

		// If we're nearing a = 1, be careful
		// Estimated step size in a is H * stepsize
		stepsize = integrator.getstepsize();
		if (data[0] + 2.0 * data[3] * stepsize > 1.0) {
			// Reduce the stepsize
			// Calculate the exact amount that the stepsize will need to be in order to get to a = 1
			// in the linear approximation
			double temp = (1.0 - data[0]) / 2.0 / data[3];
			integrator.setstepsize(0.9 * temp);
			// Eventually we'll run into the minimum step size and we'll cross the finishline
			// Also note that Hdot is usually negative, so it will typically take a little bit more than temp
			// to cross the finish line.
			// This tends to take about 20 steps to hit a = 1, which should be good enough to have
			// derivatives near a = 1 under control.
		}
	}

	// Take a look at the consistency of the final state
	check.checkfinal(data, time, params, output, status);

	return result;
}

int intfunc(double t, const double data[], double derivs[], void *params) {
	// This routine calculates the derivatives for a, phi, \dot{\phi} and \dot{H}
	// It is called by the integration routine.

	// Extract parameters
	IntParams myParams = *(IntParams *) params;

	// Call the derivatives routine in the model to calculate the derivatives appropriately
	// Note that they don't depend on time
	return myParams.getmodel().derivatives(data, derivs, myParams.getparams());

}

string getfilename(const std::string &dir, const std::string &filebase, const int padding) {
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
