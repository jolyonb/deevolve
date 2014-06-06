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
// Also need to include MCMC specific header
#include "main_MCMC.h"

using std::cout;
using std::endl;
using std::setprecision;
using std::scientific;

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
	 
    // Set up the output class. Print2Memory is used for MCMC
    Print2Memory myOutput;

    // Check that output is a go
    if (!myOutput.filesready()) {
        // End gracefully if not
        cout << "Unable to open file for output." << endl;
        return -1;
    }

    // Start timing!
    boost::timer::cpu_timer myTimer;

	// Get the priors
	ifstream priorsin;
	priorsin.open(inifile.getiniString("priorsfile", "priors.txt", "MCMC").c_str());
	vector<PARAMPRIORS> spriors;
	if(priorsin){
		while(!priorsin.eof()){
			PARAMPRIORS temp;
			priorsin >> temp.section >> temp.name >> temp.lower >> temp.upper >> temp.sigma;
			spriors.push_back(temp);
		}
		priorsin.close();
	}
	
	// Number of parameters
	int numparams = spriors.size();
	
	// Report priors info to screen
	cout << "The priors are" << endl;
	for(int n = 0; n < numparams; n++)
		cout << spriors[n].section << " " << spriors[n].name
			 << " :: " << spriors[n].lower << "\t" << spriors[n].upper << "\t" << spriors[n].sigma << endl;
	
	// See random number generator
	RNGtool.seed(time(NULL));
	
	// "temporary" holding values of the parameters
	// Values of the parameters
	double *parameters = new double[numparams];
	// Current values of the parameters
	double *current = new double[numparams];
	// Proposed values of the parameters
	double *proposed = new double[numparams];
	// Array to hold the prior info
	double *priors = new double[3 * numparams];

	// Names of the parameters
	string *names = new string[numparams];
	// Sections that the parameters are in
	string *sections = new string[numparams];	
	// Number of MCMC steps to burn
	int MCMCburninsteps = inifile.getiniInt("MCMCburninsteps", 1000, "MCMC");
	// Frequency to dump chain info to file
	int MCMCchaindumpfreq = inifile.getiniInt("MCMCchaindumpfreq", 100, "MCMC");
	// Number of MCMC steps to take
	int MCMCnumsteps = inifile.getiniInt("MCMCnumsteps", 40000, "MCMC");
	// Number of MCMC chains
	int MCMCnumchains = inifile.getiniInt("numchains", 5, "MCMC");
	// Define some useful ints and doubles
	int MCMCchainID, MCMCstep, MCMCaccept_counter;
	double lower, upper;
	
	// Populate prior array from input prior struct
	for(int n = 0; n < numparams; n++){
		names[n] = spriors[n].name;
		sections[n] = spriors[n].section;
		// priors contains the lower, upper, and sigma values
		priors[n] = spriors[n].lower;
		priors[n + numparams] = spriors[n].upper;
		priors[n + 2 * numparams] = spriors[n].sigma;
	}
	
	// Run the chains
	for(int chain = 0; chain < MCMCnumchains; chain++){
		
		// Report chain number to screen
		cout << "chain # " << chain << endl;
		// Zero the MCMC step number
		MCMCstep = 0;
		// Zero the MCMC acceptance counter
		MCMCaccept_counter = 0;
		// ID of this chain
		MCMCchainID = chain+1E4;
		
		// Start off at a random position in parameter space
		for(int param = 0; param < numparams; param++){
			lower = priors[param];
			upper = priors[param + numparams];
			parameters[param] = lower + UnitRand() * (upper - lower);
		}
	
		// Open up file to dump chain info
		ofstream MCMCchainfile;
		string filename = inifile.getiniString("chainddir", "chains", "MCMC") +"/"
						  + inifile.getiniString("chaindir", "run1", "MCMC") +"/"				
						  + inifile.getiniString("chainfileprefix", "chains", "MCMC")  
						  + "_" + Int2String(MCMCchainID) + ".dat";
		MCMCchainfile.open(filename.c_str());
		
		// Variables for holding the current and proposed likelihoods,
		// and their ratio
	    double L_current, L_proposed, LikelihoodRatio;
     	
		// Start the sampling
		while(true){
		
			// Dump the values of "parameters" into "current"
			// NOTE: for step = 0, these are the random initial points picked
			memcpy(current, parameters, numparams * sizeof(double));
			// Current value of the likelihood function
			L_current = ComputeLikelihood(inifile, names, sections, current, numparams);
			// Get some proposed parameters
			GetProposedParameters(priors, current, proposed, numparams);
			// Get the value of the likelihood with these proposed parameters		
			L_proposed = ComputeLikelihood(inifile, names, sections, proposed, numparams);
			// Compute likelihood ratio	
			LikelihoodRatio = L_proposed / L_current;
			// Decide whether to accept the proposed parameters.
			if( L_proposed >= L_current || UnitRand() < LikelihoodRatio ){
				// Increment acceptance counter
				MCMCaccept_counter++;
				// Store proposed parameters
				memcpy(parameters, proposed, numparams*sizeof(double));
			}
		
			// Dump to file after burn-in
			if(MCMCstep > MCMCburninsteps){
			
				// Only dump chain info to file every "chaindumpfreq" MCMCsteps.				
				if(MCMCstep % MCMCchaindumpfreq == 0 ){
					for(int n = 0; n < numparams; n++)
						MCMCchainfile << parameters[n] << "\t";
					MCMCchainfile << "\t" << L_current << endl;
				}
			
			}
			else{
			    // Only start the acceptance counter after the burn-in period has ended
				MCMCaccept_counter = 0;
			}
		
			// Say we've taken one more step
			MCMCstep++;
			// Halt if we've taken enough steps
			if(MCMCstep >= MCMCnumsteps) break;
			
		} // END samping while(){}		   
	
		// Write final counts for this chain
		MCMCchainfile << "# Number of samples (after burn-in): " << MCMCstep - MCMCburninsteps
		               << ", Number of acceptances: " << MCMCaccept_counter << endl;
		MCMCchainfile.close();
		cout << setprecision(4) << "Acceptance rate = " 
			 << 100 * MCMCaccept_counter / (float) (MCMCstep - MCMCburninsteps) << "%" << endl;

	} // END chain-loop

    //**********//
    // Clean up //
    //**********//

	delete parameters;
	delete current;
	delete proposed;
	delete priors;

    // Stop timing
    myTimer.stop();
    double ms = myTimer.elapsed().wall / 1e6;
    if (ms < 1e3)
        cout << setprecision(4) << "chains complete in " << ms << " milliseconds.";
    else
        cout << setprecision(4) << "chains complete in " << ms / 1000.0 << " seconds.";
    cout << endl;

	// Exit gracefully
	return result;
	
}


double ComputeLikelihood(IniReader& inifile, string *names, string *sections, double *parameters, int numparams){
	
	// Set values of the parameters (note: need to get name & section for inifile)
	for(int n = 0; n < numparams; n++)	
		inifile.setparam(names[n], sections[n], parameters[n]);
	
    // Set up the cosmological parameters 
    Parameters myParams(inifile);
	// Setup the Print2Memory class
    Print2Memory myOutput;
	// Do the evolution, and return the likelihood for the data combination
	// defined by "combinationchi"
    if ( doEvolution(inifile, myParams, myOutput, true) == 0) 
		return exp( - 0.5 * myOutput.getvalue("combinationchi", 0.0));
	else
		return 0;	
	
} // END ComputeLikelihood()
 
 
void GetProposedParameters(double *priors, double *current, double *proposed, int numparams){
	
	// Get the upper and lower bounds, and sigma on the prior for this parameter
	double upper, lower, sigma;
	// Holding variable of the proposed parameter value
	double prop;
	// Holding variable of the current parameter value
	double thisval;
	// Loop over all parameters
	for(int param = 0; param <  numparams; param++){
		
		thisval = current[param];
		lower = priors[param];
		upper = priors[param + numparams];
		sigma = priors[param + 2 * numparams];

		// This process may choose parameter values which live outside the prior range,
		// and so must repeat until a parameter is found which is inside prior range.
		while( true ){
			
			// Use Box-Muller transform to get N(0,1) -- Normally distributed number.
			prop = thisval + BoxMuller() * sigma;
			
			// Make sure the proposal is inside the prior range before getting out
			if(prop <= upper && prop >= lower) break;
			
		}
		// Place the new proposal into the array
		proposed[param] = prop;
	}
	
} // END GetProposedParameters() 