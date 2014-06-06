/*
 * main.cpp
 *
 * This is the program entry point for the program.
 *
 * The purpose of this software is to evolve scalar field dark energy models through cosmological time.
 *
 * This wrapper performs a MCMC exploration of parameter space.
 *
 * This software requires the GSL libraries and the C++ BOOST libraries (see www.boost.org)
 * These are most easily installed using a package manager (libboost-all-dev on ubuntu)
 *
 * Jolyon K. Bloomfield and Jonathan A. Pearson, March 2014
 *
 */

#include "main.h"

using namespace std;

// Structure for storing results in a vector
typedef struct expresults {
    double data[7];
} expresults;
const int numchisqds = 7;

// Random number generation tools
boost::random::mt19937 RNGtool;
boost::random::uniform_real_distribution<> RNGunit(0.0, 1.0);
boost::random::normal_distribution<> RNGnormal(0.0, 1.0);
double UnitRand(){ return RNGunit(RNGtool); }
double NormalRand(){ return RNGnormal(RNGtool); }

string Int2String(int Number) { return static_cast<ostringstream*>( &(ostringstream() << Number) )->str(); }

// A routine that prints a progressbar to screen
// Make sure to call std::cout << endl; to clear the bar after finishing
void updateprogress(float progress) {
    int barWidth = 70;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

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

	bool showprogress = inifile.getiniBool("progress", true, "MCMC");

    //**************//
    // Output class //
    //**************//

    // Set up the filenames to output
    string outputdir = inifile.getiniString("logdir", "logs", "Function");
    string basename = inifile.getiniString("runname", "run", "Function");

    // Go and find our appropriate file name (using 4 digit numbers as the default)
    string outputname = getfilename(outputdir, basename, "", inifile.getiniInt("numberpad", 4, "Function"), true) + ".dat";

    // Set up the output class. Print2Memory is used for MCMC
    Print2Memory myOutput;

    // Check that output is a go
    std::ofstream outputstream(outputname.c_str());
    if (!outputstream.is_open()) {
        // End gracefully if not
        cout << "Unable to open file for output." << endl;
        return -1;
    }


    //*****************//
    // Set up the MCMC //
    //*****************//

    int numsteps = inifile.getiniInt("numsteps", 10000, "MCMC");
    int burnsteps = inifile.getiniInt("burnsteps", 1000, "MCMC");
    string priorfile = inifile.getiniString("priorfile", "mcmc.txt", "MCMC");



    // Report what we're doing to screen


    cout << "Outputting to " << outputname << endl;



    // Progressbar stuff
    float progress = 0.0;
    int barcount = 0;

    //*******************//
    // Perform the sweep //
    //*******************//

    // Start timing!
    boost::timer::cpu_timer myTimer;

    // Loop through the parameter space
    // Second parameter first (only runs through once if not being used)
    for (int i2 = 0; i2 < numsteps2; i2++) {
        // If we are doing a two parameter scan, update the second parameter here
        stepper2 = lower2 + stepsize2 * i2;
        if (numparams == 2)
            inifile.setparam(param2, section2, stepper2);

        // First parameter next
        for (int i1 = 0; i1 < numsteps1; i1++) {
            // Update the first parameter here
            stepper1 = lower1 + stepsize1 * i1;
            inifile.setparam(param1, section1, stepper1);

            // Set up the cosmological parameters (done here in case something significant changed in the sweeping)
            Parameters myParams(inifile);

            // Do the evolution
            result = doEvolution(inifile, myParams, myOutput, true);

            // Interpret the result of the evolution
            if (result == 0) {
                // Everything was successful. Now we can save the results!

                // Store the parameter values
                parameter1.push_back(stepper1);
                if (numparams == 2) parameter2.push_back(stepper2);

                // Extract the chi^2 values
                filling.data[0] = myOutput.getvalue("WMAPchi", 0.0);
                filling.data[1] = myOutput.getvalue("PLANCKchi", 0.0);
                filling.data[2] = myOutput.getvalue("SNchi", 0.0);
                filling.data[3] = myOutput.getvalue("Hubblechi", 0.0);
                filling.data[4] = myOutput.getvalue("BAOtotalchi", 0.0);
                filling.data[5] = myOutput.getvalue("BAOtotalrchi", 0.0);
                filling.data[6] = myOutput.getvalue("combinationchi", 0.0);
                // Plop that on the stack too!
                chisquareds.push_back(filling);

            }

            // Update the progress bar every 20 cycles
            if (showprogress && ++barcount >= 20) {
                barcount = 0;
                progress = (float) (i2 * numsteps1 + i1 + 1) / (float) totnumsteps;
                updateprogress(progress);
            }

        }

        // Update the progress bar (to make sure it shows 100% at the end)
        barcount = 0;
        progress = (float) ((i2 + 1) * numsteps1) / (float) totnumsteps;
        if (showprogress) updateprogress(progress);
    }
    // Clear the progressbar
    if (showprogress) std::cout << std::endl;


    //**********//
    // Clean up //
    //**********//

    outputstream.close();

    // Stop timing
    myTimer.stop();
    double ms = myTimer.elapsed().wall / 1e6;
    if (ms < 1e3)
        cout << setprecision(4) << "MCMC complete in " << ms << " milliseconds.";
    else
        cout << setprecision(4) << "MCMC complete in " << ms / 1000.0 << " seconds.";
    cout << endl;

	// Exit gracefully
	return 0;
}
