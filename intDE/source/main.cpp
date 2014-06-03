/*
 * main.cpp
 *
 * This is the program entry point for the program.
 *
 * It acts as a wrapper around the evolution routines. All that this wrapper does is to set
 * up the appropriate input and output objects.
 *
 * The purpose of this software is to evolve scalar field models through cosmological time.
 *
 * This software requires the GSL libraries.
 *
 * This software requires the C++ BOOST libraries (see www.boost.org)
 * These are most easily installed using a package manager (libboost-all-dev on ubuntu)
 *
 * Jolyon K. Bloomfield and Jonathan A. Pearson, March 2014
 *
 */

#include "main.h"

using std::cout;
using std::endl;
using std::setprecision;

// Our program entry point
// This entry point is just a wrapper around the evolution routines.
// It sets up the input parameters as well as the output file, and otherwise just calls the routines.
int main(int argc, char* argv[]) {

    int result; // Just a number for return values

    //******************//
    // Input parameters //
    //******************//

	// Read the input parameters file, named "params.ini" by default.
	// If there is a command line argument, assume that it is the filename for the parameters file.
	IniReader inifile;
	if (argc > 1)
		inifile.read(argv[1]);
	else
		inifile.read("params.ini");

	// Values in the ini file can be set using the following command:
    // inifile.setparam("desiredh", "Cosmology", 1);
	// First parameter is the key name, the second is the section name, and the third is the value, either an integer, string or double

	// Check to see if we're doing an individual evolution or sweeping over a parameter
	if (inifile.getiniBool("sweeping", false, "Function") == true) {
	    result = doSweep(inifile);
	} else {
	    result = doSingleEvolution(inifile);
	}

	// Exit gracefully
	return result;
}

string getfilename(const std::string &dir, const std::string &filebase, const std::string &postbase, const int padding) {
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
		std::ostringstream convert;
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
		if (exists(dir + "/" + filebase + filenum + postbase + ".dat"))
			continue;

		// If we got to here, we have a unique filename; return it
		return dir + "/" + filebase + filenum;
	}

	// We really shouldn't get here, but parsers like making sure there's a return
	return "error";

}

// Routine to perform a single evolution
int doSingleEvolution(IniReader &inifile) {

    int result;

    // Set up the cosmological parameters
    Parameters *myParams = new Parameters(inifile);

    //**************//
    // Output class //
    //**************//

    // Set up the filenames to output
    string outputdir = inifile.getiniString("logdir", "logs", "Function");
    string basename = inifile.getiniString("runname", "run", "Function");
    string postname = inifile.getiniString("postname", "d", "Function");

    // Go and find our appropriate file name (using 4 digit numbers as the default)
    string outputname = getfilename(outputdir, basename, postname, inifile.getiniInt("numberpad", 4, "Function"));

    // Set up the output class
    std::string parsestring = inifile.getiniString("outputclass", "BasicDump", "Function");
    Output *myOutput;
    if (parsestring == "BasicDump")
        myOutput = new BasicDump(outputname, postname);
    else if (parsestring == "Print2Memory")
        myOutput = new Print2Memory(outputname, postname);
    else
        myOutput = new BasicDump(outputname, postname);    // BasicDump is the default

    // Check that output is a go
    if (!myOutput->filesready()) {
        // End gracefully if not
        cout << "Unable to open files for output." << endl;
        delete myOutput;
        delete myParams;
        return -1;
    }

    // Initialize vectors to store Hubble and redshift data
    vector<double> hubble;
    vector<double> redshift;


    //*******************//
    // Do the evolution! //
    //*******************//

    // Print some stuff to the screen
    cout << "Beginning evolution." << endl;
    cout << "Outputting to " << outputname << endl;

    // Start timing!
    boost::timer::cpu_timer myTimer;

    // Do the evolution
    result = doEvolution(inifile, *myParams, *myOutput, redshift, hubble);

    // Stop timing
    myTimer.stop();

    // Interpret the result of the evolution
    if (result == 0) {
        // Success!

        // Print a nice message
        myOutput->printfinish(myTimer.elapsed().wall / 1e6);
        cout << setprecision(4) << "Evolution complete in " << myTimer.elapsed().wall / 1e6 << " milliseconds." << endl;

        // Perform postprocessing if specified in the options
        if (inifile.getiniBool("postprocess", false, "Function") == true) {
            // Start timing!
            boost::timer::cpu_timer myPostTimer;
            cout << "Beginning postprocessing." << endl;

            // Perform postprocessing
            result = PostProcessing(inifile, *myParams, *myOutput, redshift, hubble);

            // Stop timing
            myPostTimer.stop();

            // Check for result of postprocessing
            if (result == 1) {
                cout << "Error integrating distance measures; terminating." << endl;
            }
            else if (result == 0) {
                // Report done
                cout << "Postprocessing completed in " << setprecision(4) << myPostTimer.elapsed().wall / 1e6 << " milliseconds." << endl;
            }
        }
    }
    else if (result == -1) {
        // Initialization error
        cout << "Error initializing model; terminating." << endl;
    }
    else if (result == 1) {
        // Integration error
        cout << "Integration error; terminating." << endl;
    }
    else if (result == 2) {
        // NAN error
        cout << "NAN error; terminating." << endl;
    }
    else if (result == 3) {
        // Did not get t a = 1 error
        cout << "Did not evolve to a = 1 within appropriate conformal time; terminating." << endl;
    }
    else if (result == 4) {
        // Invalid state error
        cout << "Invalid state reached; terminating." << endl;
    }


    //**********//
    // Clean up //
    //**********//

    // Do the final print on the output class
    myOutput->printfinal("modelOmegaR");

    // No memory leaks!
    delete myOutput;
    delete myParams;

    return 0;

}

// Routine to sweep over a parameter
int doSweep(IniReader &inifile) {

    int result;

    //**************//
    // Output class //
    //**************//

    // Set up the filenames to output
    string outputdir = inifile.getiniString("logdir", "logs", "Function");
    string basename = inifile.getiniString("runname", "run", "Function");
    string postname = inifile.getiniString("postname", "d", "Function");

    // Go and find our appropriate file name (using 4 digit numbers as the default)
    string outputname = getfilename(outputdir, basename, "", inifile.getiniInt("numberpad", 4, "Function"));
    string likelihood = outputname + postname + ".dat";

    // Set up the output class -- Print2Memory is used for sweeps
    Print2Memory *myOutput = new Print2Memory(outputname, "");

    // Check that output is a go
    if (!myOutput->filesready()) {
        // End gracefully if not
        cout << "Unable to open file for output." << endl;
        delete myOutput;
        return -1;
    }

    // Initialize vectors to store Hubble and redshift data
    vector<double> hubble;
    vector<double> redshift;


    //******************//
    // Set up the sweep //
    //******************//
    string param = inifile.getiniString("param", "phi0", "Sweep");
    string section = inifile.getiniString("section", "Cosmology", "Sweep");
    double lower = inifile.getiniDouble("lower", -1.0, "Sweep");
    double upper = inifile.getiniDouble("upper", 1.0, "Sweep");
    int numsteps = inifile.getiniDouble("steps", 200, "Sweep");
    double stepsize = (upper - lower) / static_cast<double> (numsteps);

    // Storage for various values
    vector<double> parameter;
    vector<expresults> chisquareds;
    expresults filling;
    // Allocate storage space
    parameter.reserve(numsteps);
    chisquareds.reserve(numsteps);


    //*******************//
    // Perform the sweep //
    //*******************//

    // Print some stuff to the screen
    cout << "Sweeping over " << param << " from " << lower << " to " << upper << " in " << numsteps << " steps." << endl;
    cout << "Outputting to " << outputname << endl;

    // Start timing!
    boost::timer::cpu_timer myTimer;

    // Loop through the parameterspace
    for (double stepper = lower; stepper <= upper; stepper += stepsize) {

        // Set the parameters in the inireader
        inifile.setparam(param, section, stepper);

        // Set up the cosmological parameters (done here in case something significant changed)
        Parameters *myParams = new Parameters(inifile);

        // Do the evolution
        result = doEvolution(inifile, *myParams, *myOutput, redshift, hubble);

        // Interpret the result of the evolution
        if (result == 0) {
            // Perform postprocessing
            result = PostProcessing(inifile, *myParams, *myOutput, redshift, hubble);
            if (result == 0) {
                // Everything was successful. Now we can save the results!
                // Add the parameter value
                parameter.push_back(stepper);
                // Populate the filling structure
                filling.data[0] = myOutput->getvalue("WMAPchi", -1.0);
                filling.data[1] = myOutput->getvalue("PLANCKchi", -1.0);
                filling.data[2] = myOutput->getvalue("SNchi", -1.0);
                filling.data[3] = myOutput->getvalue("Hubblechi", -1.0);
                filling.data[4] = myOutput->getvalue("6dFGSchi", -1.0);
                filling.data[5] = myOutput->getvalue("SDSSchi", -1.0);
                filling.data[6] = myOutput->getvalue("SDSSRchi", -1.0);
                filling.data[7] = myOutput->getvalue("WiggleZchi", -1.0);
                filling.data[8] = myOutput->getvalue("BOSSDR9chi", -1.0);
                filling.data[9] = myOutput->getvalue("BOSSDR11chi", -1.0);
                // Combine data sets: WMAP, SN, SDSSR, WiggleZ, BOSSDR9
                filling.data[10] = filling.data[0] + filling.data[2] + filling.data[6] + filling.data[7] + filling.data[8];
                // Plop that on the stack too!
                chisquareds.push_back(filling);

                // Print the chi^2 values to file, as well as the parameter
                myOutput->printfinal(param);
            }
        }

        // Clean up
        delete myParams;

    }


    //****************//
    // Postprocessing //
    //****************//

    // We have the chi^2 values for all experiments for each value of the parameter that we've scanned over
    // What we want to do now is to convert these chi^2 values into likelihood values
    // We also do a combined likelihood, throwing different pieces together

    // Storage
    vector<expresults> likelihoods;
    numsteps = chisquareds.size(); // Just in case the computation from above was wrong?
    likelihoods.reserve(numsteps);

    // Iterate over everything, calculating the likelihood L = e^{-chi2/2}
    for (int i = 0; i < numsteps; i++) {
        filling = chisquareds[i];
        for (int j = 0; j < 11; j++) {
            filling.data[j] = exp(- 0.5 * filling.data[j]);
        }
        likelihoods.push_back(filling);
    }

    // Find the highest likelihood values, so that we can normalise using them
    for (int j = 0; j < 11; j++) {
        filling.data[j] = 0;
    }
    for (int i = 0; i < numsteps; i++) {
        for (int j = 0; j < 11; j++) {
            if (likelihoods[i].data[j] > filling.data[j])
                filling.data[j] = likelihoods[i].data[j];
        }
    }

    // Perform the normalization
    for (int i = 0; i < numsteps; i++)
        for (int j = 0; j < 11; j++)
            likelihoods[i].data[j] /= filling.data[j];

    // Output the likelihood data to file!
    std::ofstream outputstream(likelihood.c_str());
    outputstream << std::scientific << setprecision(8) << "# " << param << "\tWMAP\tPLANCK\tSN\tHubble\t6dFGS\tSDSS\tSDSSR\tWiggleZ\tBOSSDR9\tBOSSDR11\tCombined" << endl;
    for (int i = 0; i < numsteps; i++) {
        outputstream << parameter[i];
        for (int j = 0; j < 11; j++)
            outputstream << "\t" << likelihoods[i].data[j];
        outputstream << endl;
    }
    outputstream.close();

    //**********//
    // Clean up //
    //**********//

    // Stop timing
    myTimer.stop();
    cout << setprecision(4) << "Sweep complete in " << myTimer.elapsed().wall / 1e6 << " milliseconds." << endl;

    // No memory leaks!
    delete myOutput;

    return 0;

}
