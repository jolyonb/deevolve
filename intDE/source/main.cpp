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

// All the includes are hidden in the header file
#include "main.h"

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

	// Now that we have the ini file read in, we can change values programmatically as follows.
	// Start by grabbing the tree data
	// boost::property_tree::ptree initdata = inifile.getdata();
	// Change a value
	// initdata.put("Cosmology.Hubbleh", "0.7");
	// And send it back!
	// inifile.setdata(initdata);
	// How easy is that?

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
	double H0;
	result = doEvolution(inifile, *myParams, *myOutput, redshift, hubble, H0);

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

            // Lodge a call to the parameters class asking it to update it's information
            myParams->updateinfo(H0);

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

    delete myOutput;
    delete myParams;

	// Exit gracefully
	return 0;
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
		if (exists(dir + "/" + filebase + filenum + postbase + ".dat"))
			continue;

		// If we got to here, we have a unique filename; return it
		return dir + "/" + filebase + filenum;
	}

	// We really shouldn't get here, but parsers like making sure there's a return
	return "error";

}
