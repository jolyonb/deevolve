#include "main.h"

struct DATA{
	double x,y;
};
struct PFIT{
	double m;
	double c;
};

struct PARAMP{
	string name;
	double lower;
	double upper;
	double sigma;
};

struct RPS{
	int chainID;
	string chaindir;
	string chainfileprefix;
	int numMCMCsteps;
	int burninsteps;
	int numchains;
	int chaindumpfreq;
};


double ComputeLikelihood(double *params, vector<DATA> &data);
void GetProposedParameters(double *priors, double *current, double *proposed, int numparams);
void runchain(struct RPS &runparams, vector<DATA> &data, vector<PARAMP> priors);

int main(){

	// Initialize the run parameters
	RPS runparams;
	
	// How many MCMC chains do we want?
	runparams.numchains = 5;
	// Output directory of the chains
	runparams.chaindir = "chains/";
	// File name-prefix of this chain
	runparams.chainfileprefix = "chain";
	// Number of steps to run MCMC for
	runparams.numMCMCsteps = 40000;
	// Number of steps which are to be burnt
	runparams.burninsteps = 1000;
	// After how many MCMC steps should chain info be dumped to file?
	runparams.chaindumpfreq = 10;
	// File name of the test data file
	string datafile = "testdata.dat";
	// File name of priors file
	string priorsfile = "priors.txt";

	//////////////////////////////////
	// Some stuff specific to this MCMC
	//////////////////////////////////
	
	// Get the data we want to fit
	ifstream indata;
	indata.open(datafile.c_str());
	vector<DATA> data;
	if(indata){
		while(!indata.eof()){
			DATA temp;
			indata >> temp.x >> temp.y;
			data.push_back(temp);
		}
		indata.close();
	}
	
	// Get the priors
	ifstream gps;
	gps.open(priorsfile.c_str());
	vector<PARAMP> priors;
	if(gps){
		while(!gps.eof()){
			PARAMP temp;
			gps >> temp.name >> temp.lower >> temp.upper >> temp.sigma;
			priors.push_back(temp);
		}
		gps.close();
	}
	cout << "The priors are" << endl;
	for(int n = 0; n < (int) priors.size(); n++)
		cout << priors[n].name << " :: " << priors[n].lower << "\t" << priors[n].upper << "\t" << priors[n].sigma << endl;
	
	// Seed the random number generator
	// srand (time(NULL));
	RNGtool.seed(time(NULL));
	// Note: if you're running multiple chains in separate processes starting around the same time, make sure to seed them using
	// different seeds!
	
	// Run the chains!
	for(int chain = 0; chain < runparams.numchains; chain++){
		runparams.chainID = 1E4 + chain;
		runchain(runparams, data, priors);
	}
	
	cout << "Done" << endl;

	return 0;
	
} // END main();

void runchain(struct RPS &runparams, vector<DATA> &data, vector<PARAMP> spriors){
	
	// Zero the MCMC step number
	int MCMCstep = 0;
	// Zero the acceptance counter
	int MCMCaccept_counter = 0;
	// Inherit the number of MCMC steps
	int MCMCnumsteps = runparams.numMCMCsteps;
	// Inherit the fraction of MCMC steps to burn.
	int MCMCburninsteps = runparams.burninsteps;
	// How often to keep samples?
	int MCMCchaindumpfreq = runparams.chaindumpfreq;
	// Output stream for chaininfo
	ofstream MCMCchainfile;
	// The likelihood variables
	double L_current, L_proposed, LikelihoodRatio;
	
	// Get the number of parameters
	int numparams = spriors.size();
	// Variable holding the values of the parameters	
	double *parameters = new double[numparams];
	// Create array to hold the current values of the parameters
	double *current = new double [numparams];
	// Create array to hold the proposed values of the parameters
	double *proposed = new double [numparams];	
	// Array holding prior info
	double *priors = new double[3 * numparams];

	// Populate prior array from input prior struct
	for(int n = 0; n < numparams; n++){
		priors[n] = spriors[n].lower;
		priors[n + numparams] = spriors[n].upper;
		priors[n + 2 * numparams] = spriors[n].sigma;
	}
	
	// Start off at a random position in parameter space
	double lower, upper;
	for(int param = 0; param < numparams; param++){
		lower = priors[param];
		upper = priors[param + numparams];
		parameters[param] = lower + UnitRand() * (upper - lower);
	}
	
	// Open up file to dump chain info
	string filename = runparams.chaindir + runparams.chainfileprefix + "_" + Int2String(runparams.chainID) + ".dat";
	MCMCchainfile.open(filename.c_str());
	
	// Start the sampling
	while(true){
		
		// Dump the values of "parameters" into "current"
		// NOTE: for step = 0, these are the random initial points picked
		memcpy(current, parameters, numparams * sizeof(double));
		
		// Current value of the likelihood function
		L_current = ComputeLikelihood(current, data);
		
		// Get some proposed parameters
		GetProposedParameters(priors, current, proposed, numparams);
		
		// Get the value of the likelihood with these proposed parameters		
		L_proposed = ComputeLikelihood(proposed, data);
			
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
	}

	// Write final counts
	MCMCchainfile << "# Number of samples (after burn-in): " << MCMCstep - MCMCburninsteps
	               << ", Number of acceptances: " << MCMCaccept_counter << endl;
	MCMCchainfile.close();

	cout << setprecision(4) << "Acceptance rate = " << 100 * MCMCaccept_counter / (float) (MCMCstep - MCMCburninsteps) << "%" << endl;

	delete parameters;
	delete current;
	delete proposed;
	delete priors;
	
} // END runchain()

double ComputeLikelihood(double *params, vector<DATA> &data){
	
	double test_m = params[0];
	double test_c = params[1];
	double var = 0.0;
	double testy, diff;
	for(int n = 0; n < (int) data.size(); n++ ){
		testy = test_m * data[n].x + test_c;
		diff = data[n].y - testy;
		var += diff * diff;
	}
	
	return exp(- 0.5 * var );
	
}// END ComputeLikelihood()

void GetProposedParameters(double *priors, double *current, double *proposed, int numparams){
	
	// Get the upper and lower bounds, and sigma on the prior for this parameter
	double upper, lower, sigma;
	// temporary holding variable of the proposed parameter value
	double prop;
	// temporary holding variable of the current parameter value
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
