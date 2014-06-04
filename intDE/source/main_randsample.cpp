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
	string sampdensity = outputname + "s" + ".dat";
	
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
	string section = inifile.getiniString("section", "Cosmology", "Sweep");
	int numsteps;

	/*
	string param = inifile.getiniString("param", "phi0", "Sweep");
    double lower = inifile.getiniDouble("lower", -1.0, "Sweep");
    double upper = inifile.getiniDouble("upper", 1.0, "Sweep");
    int numsteps = inifile.getiniDouble("steps", 200, "Sweep");
    double stepsize = (upper - lower) / static_cast<double> (numsteps);
*/
	// JAP
	// Get file name containing the sweep parameters
	string sweepsfile = inifile.getiniString("sweepsfile", "sweeps.txt", "Sweep");
	// Vector to hold sweep parameters
	vector<UPARAMS> iparams;
	// Open up sweeps file
	std::ifstream UI;
	UI.open(sweepsfile);
	bool catchsingle = false;
	// Read in sweeps file; store in iparams
	if(UI){
		while(!UI.eof()){
			UPARAMS ptemp;
			UI >> ptemp.name >> ptemp.lower >> ptemp.upper >> ptemp.stepsize;
			ptemp.numsteps = floor((ptemp.upper - ptemp.lower) / ptemp.stepsize) + 1;
			if(ptemp.name.substr(0,1)!="#") {
				iparams.push_back(ptemp);
				if(ptemp.numsteps==1)catchsingle = true;
			}
		}
		UI.close();
	}
	// Report to screen
	int totnumsteps = 1;
	for(int n = 0; n < iparams.size(); n++){
		cout << "Sweeping over ";
		cout << iparams[n].name << " from " << iparams[n].lower << " to  ";
		cout << iparams[n].upper << ", with step-size " << iparams[n].stepsize ;
		cout << " (number of steps = " << iparams[n].numsteps << ")" << endl;
		totnumsteps*=iparams[n].numsteps;
	}
	cout << "Total number of samples = " << totnumsteps << endl;
	cout << "Outputting to " << outputname << endl;
	// !JAP
	
    // Storage for various values
    vector<double> parameter1;
    vector<double> parameter2;	
    vector<expresults> chisquareds;
    expresults filling;
    // Allocate storage space
    parameter1.reserve(totnumsteps);
    parameter2.reserve(totnumsteps);	
    chisquareds.reserve(totnumsteps);


    //*******************//
    // Perform the sweep //
    //*******************//


    // Start timing!
    boost::timer::cpu_timer myTimer;

    // Loop over the desired number of samples
	double paramval[2];
	paramval[0] = iparams[0].lower + rand()/(double)RAND_MAX * (iparams[0].upper - iparams[0].lower);
	paramval[1] = iparams[1].lower + rand()/(double)RAND_MAX * (iparams[1].upper - iparams[1].lower);
	double combinedlike;
	double maxlike = 0.0;
	double acceptthresh = inifile.getiniDouble("acceptthresh", 0.1, "Sweep");
	int numsamples = inifile.getiniInt("numsamples", 10, "Sweep");
	double burnfrac = inifile.getiniDouble("burnfrac", 0.1, "Sweep");
	
	int burnsamples = int( burnfrac * numsamples );
	int shakeup = inifile.getiniInt("shakeup", 10, "Sweep");
	bool getnewvalues[2] = {true,true};

	for (int sample = 0; sample < numsamples; sample++) {
        
		// Set the parameters in the inireader
		
		inifile.setparam(iparams[0].name, section, paramval[0]);
        inifile.setparam(iparams[1].name, section, paramval[1]);

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
                parameter1.push_back(paramval[0]);
				parameter2.push_back(paramval[1]);
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
                myOutput->printfinal(iparams[0].name);
                myOutput->printfinal(iparams[1].name);				
				
				// Get the combined likelihood for this parameter combination
				combinedlike = exp(- 0.5 * filling.data[10]);
				
				// We only ever want to "not" modify our parameters after a burn-in period
				if( sample > burnsamples ){
					// If the current parameter choice yielded a smaller chi2 than the previous choice,
					// first remember this chi2, then say that we dont want new parameter values.
					if( combinedlike >= maxlike ){
						maxlike = combinedlike;
						for(int i = 0; i < 2; i++){
							if( rand()/(double)RAND_MAX < acceptthresh ) {
								getnewvalues[i] = true;
							}
							else{
								getnewvalues[i] = false;
							}
						}
					}
					else{
						for(int i = 0; i < 2; i++){
							getnewvalues[i] = true;
						}
					}
					// Force change of parameters every "shakeup" samples
					if( sample % shakeup == 0 ){
						for(int i = 0; i < 2; i++){
							getnewvalues[i] = true;
						}
					}
				} // END if(sample > burnsamples){}
				
				
				
				// If after all this we want new values, get them randomly here.
				for(int i = 0; i < 2; i++){
					if(getnewvalues[i]){
						// THis should probably be modified: how far to move MUST be 
						// dependnent on the likelihood...?
						paramval[i] = iparams[i].lower + rand()/(double)RAND_MAX * (iparams[i].upper - iparams[i].lower);
					}
				}
				 
				
            } // END if (result == 0){}
        } // END if (result == 0){}
		
		

        // Clean up
        delete myParams;
	} // END sample


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
    outputstream << std::scientific << setprecision(8) << "# " ;
	outputstream << iparams[0].name << " " << iparams[1].name << " ";
	outputstream << "\tWMAP\tPLANCK\tSN\tHubble\t6dFGS\tSDSS\tSDSSR\tWiggleZ\tBOSSDR9\tBOSSDR11\tCombined" << endl;
	double prevparam1 = parameter1[0];
    for (int i = 0; i < numsteps; i++) {
		if(parameter1[i]!=prevparam1 && !catchsingle) {outputstream << endl; prevparam1 = parameter1[i];}
        outputstream << parameter1[i] << " " << parameter2[i];
        for (int j = 0; j < 11; j++)
            outputstream << "\t" << likelihoods[i].data[j];
        outputstream << endl;
    }
    outputstream.close();


	// Create histogram of the combined likelihood
	int nbins[2];
	nbins[0] = inifile.getiniInt("nbins", 10, "Sweep");;
	nbins[1] = inifile.getiniInt("nbins", 10, "Sweep");;
	double dp[2];
	for(int n=0; n < 2; n++){
		dp[n] = (iparams[n].upper - iparams[n].lower)/nbins[n];
	}
	double p1_min, p1_max;
	double p2_min, p2_max;
	
	vector<vector <int> > bin;
	
	for(int i = 0; i < nbins[0]; i++){
		vector<int>dum;
		for(int j = 0; j < nbins[1]; j++){
			dum.push_back(0);
		}
		bin.push_back(dum);
	}
	double dummycounter = 0;
	for(int i = 0; i < nbins[0]; i++){
		// lower bound on bin
		p1_min = iparams[0].lower + i*dp[0];
		p1_max = p1_min + dp[0];
		for(int j = 0; j < nbins[1]; j++){
			p2_min = iparams[1].lower + j*dp[1];
			p2_max = p2_min + dp[1];
			for(int n=0; n < numsteps; n++){
				if( p1_min<parameter1[n] && parameter1[n]<p1_max ){
					if(p2_min<parameter2[n] && parameter2[n]<p2_max){
						bin[i][j]++;
						dummycounter ++;
					}
				}
			}
		}
	}
	
    std::ofstream sdout(sampdensity.c_str());
	for(int i = 0; i < nbins[0]; i++){
		for(int j = 0; j < nbins[1]; j++){
			sdout << iparams[0].lower + i*dp[0] << " ";
			sdout << iparams[1].lower + j*dp[1] << " ";
			sdout << bin[i][j] << endl;
		}
		sdout << endl;
	}
	sdout.close();

    //**********//
    // Clean up //
    //**********//

    // Stop timing
    myTimer.stop();
    cout << setprecision(4) << "Sweep complete in " << myTimer.elapsed().wall / 1e6 << " milliseconds.";
	cout << endl;

    // No memory leaks!
    delete myOutput;

    return 0;

}
