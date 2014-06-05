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

// Routine to perform a sweep over a parameter
int doSweep(IniReader &inifile);

// Structure for storing results from experiments
typedef struct expresults {
    double data[11];
} expresults;

// JAP
struct UPARAMS{
    string name;
    double lower;
    double upper;
    double stepsize;
    int numsteps;
};

// !JAP

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
    double stepsize = inifile.getiniDouble("stepsize", 0.01, "Sweep");
    int numsteps = floor((upper - lower) / stepsize) + 1;

	// JAP
	// Get file name containing the sweep parameters
	string sweepsfile = inifile.getiniString("sweepsfile", "sweeps.txt", "Sweep");
	// Vector to hold sweep parameters
	vector<UPARAMS> iparams;
	vector<SWP> dparams;
	// Open up sweeps file
	std::ifstream UI;
	UI.open(sweepsfile);
	// Read in sweeps file; store in iparams
	if(UI){
		while(!UI.eof()){
			UPARAMS ptemp;
			UI >> ptemp.name >> ptemp.lower >> ptemp.upper >> ptemp.stepsize;
			ptemp.numsteps = floor((ptemp.upper - ptemp.lower) / ptemp.stepsize) + 1;
			if(ptemp.name.substr(0,1)!="#") {
				iparams.push_back(ptemp);
			/*	for(double val = ptemp.lower; val <=ptemp.upper; val +=ptemp.stepsize){
					SWP dtemp;
					dtemp.name = ptemp.name;
					dtemp.value = val;
					dparams.push_back(dtemp);
				}
				*/
			}
		}
		UI.close();
	}
	int numparams = iparams.size();
	
	
	vector< vector<double> > pcmbs;
	if(numparams==2){
		for(int n=0; n < iparams[0].numsteps; n++){
			vector<double> thispar;
			thispar.push_back
			for(int nn=0; nn < iparams[1].numsteps; nn++){
				
			}
		}
	}
	
	/*
	// Create set of vectors;
	// Each vector contains a vector of all combinations of the values of a
	// particular parameter.
	vector< vector<double> > eachparam;
	for(int n = 0; n < numparams; n++){
		vector<double> thisparam;
		for(int m = 0; m < iparams[n].numsteps; m++){
			thisparam.push_back(iparams[n].lower+m*iparams[n].stepsize);
		}
		eachparam.push_back(thisparam);
	}
	for(int PARAMID = 0; PARAMID < numparams; PARAMID++){
		for(int PARAMVAL = 0; PARAMVAL < iparams[PARAMID].numsteps; PARAMVAL++){
				cout << eachparam[PARAMID][PARAMVAL] << " ";
		}
		cout << endl;
	}
	
	for(int PARAMID = 0; PARAMID < numparams; PARAMID++){
		for(int PI2 = PARAMID+1; PI2 < numparams; PI2++){
			for(int PARAMVAL = 0; PARAMVAL < iparams[PI2].numsteps; PARAMVAL++){

				cout << eachparam[PI2][PARAMVAL] << " ";
			}
		}
		cout << endl;
	}
	
	int numcombs = 1;
	for(int n = 0; n < numparams; n++){
		numcombs*=iparams[n].numsteps;
	}
	
	vector< vector<double> > paramcombs;
	for(int c = 0; c < numcombs; c++){
		vector<double> thiscomb;
		
		for(int n = 0; n < numparams; n++){
			for(int nn=0; nn < iparams[n].numsteps; nn++){
				thiscomb.push_back(eachparam[n][nn]);
			}
		}
		
		paramcombs.push_back(thiscomb);
	}
	
	int n = 0;
	while( n < iparams[0].numsteps ){
		int i = 0;
		for(int k = 0; k < numparams; k++ ){
			i++;
			if( i == iparams[k].numsteps ) i = 0;
		}
		n++;	
	}
	*/
	// Report to screen
	for(int n = 0; n < numparams; n++){
		cout << "Sweeping over ";
		cout << iparams[n].name << " from " << iparams[n].lower << " to ";
		cout << iparams[n].upper << ", with step-size " << iparams[n].stepsize ;
		cout << " (number of steps = " << iparams[n].numsteps << ")" << endl;
	}
	
	cout << endl;
//	cout << "dp size = " << dparams.size() <<endl;
	
	
	/*
	for(int n = 0; n < iparams[0].numsteps; n++){
		vector<SWP> Vals;
		SWP tem;
		//cout << dparams[n].value << " ";
		tem.name = dparams[n].name;
		tem.value = dparams[n].value;
		Vals.push_back(tem);

		for(int nn=n; nn < dparams.size(); nn++){
			if(dparams[n].name!=dparams[nn].name){
		//		cout << dparams[nn].value << " ";
				Vals.push_back(dparams[nn].value);
			}
		}
		cout << "(";
		for(int nn=0; nn < Vals.size(); nn++){
			cout << Vals[n] << " ";
		}
		cout << ")" << endl;
	}
	*/
	  
	// !JAP
	
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
	
	// Start timing!
	boost::timer::cpu_timer myTimer;
	
	for(int n = 0; n < numparams; n++){
		// Print some stuff to the screen
		//cout << "Doing " << iparams[n].name << " from " << iparams[n].lower << " to " << iparams[n].upper << " in " << iparams[n].numsteps << " steps." << endl;
		// Loop through the parameterspace
		for (double s = iparams[n].lower; s <= iparams[n].upper; s += iparams[n].stepsize) {
		
			    inifile.setparam(iparams[n].name, section, s);

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
			            parameter.push_back(s);
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
			            myOutput->printfinal("wnaught");
			        }
			    }
			    // Clean up
			    delete myParams;
			
			} // END s-loop
		
	} // END n-loop

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
	outputstream << std::scientific << setprecision(8) << "# ";
	for(int n = 0; n < numparams; n++)
    	outputstream << std::scientific << setprecision(8) << iparams[n].name << "\t";
	
	outputstream << std::scientific << setprecision(8) << "\tWMAP\tPLANCK\tSN\tHubble\t6dFGS\tSDSS\tSDSSR\tWiggleZ\tBOSSDR9\tBOSSDR11\tCombined" << endl;
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
