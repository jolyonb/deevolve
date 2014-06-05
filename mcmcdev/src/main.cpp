
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
};

struct RPS{
	int chainID;
	string chaindir;
	string chainfileprefix;
	int numMCMCsteps;
	double burninfrac;
	int numchains;
};

struct MarkovChain{
	bool run;
	bool accept;
	int step;
	int numsteps;
	int accept_counter;
	int burnintime;
	double burninfrac;
	ofstream chainfile;
};

string Int2String(int Number) {
    return static_cast<ostringstream*>( &(ostringstream() << Number) )->str();
}


double ComputeLikelihood(vector<double> params, vector<DATA> &data);
vector<double> GetProposedParameters(  vector<PARAMP> priors, vector<double> current, double L_current);
void runchain(struct RPS &runparams, vector<DATA> &data, vector<PARAMP> priors);

int main(){

	// Initialize the run parameters
	RPS runparams;
	
	// How many MCMC chains do we want?
	runparams.numchains = 2;
	// Output directory of the chains
	runparams.chaindir = "chains/";
	// File name-prefix of this chain
	runparams.chainfileprefix = "chain";
	// Number of steps to run MCMC for
	runparams.numMCMCsteps = 40000;
	// Fraction of total steps which are to be burnt
	runparams.burninfrac = 0.1;
	
	// File name of the test data file
	string datafile = "testdata.dat";
	// File name of priors file
	string priorsfile = "priors.txt";

	//////////////////////////////////
	// Some stuff specific to this MCMC
	//////////////////////////////////
	
	// Get the data we want to fit
	ifstream indata;
	indata.open(datafile);
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
	gps.open(priorsfile);
	vector<PARAMP> priors;
	if(gps){
		while(!gps.eof()){
			PARAMP temp;
			gps >> temp.name >> temp.lower >> temp.upper;
			priors.push_back(temp);
		}
		gps.close();
	}
	cout << "The priors are" << endl;
	for(int n = 0; n < priors.size(); n++)
		cout << priors[n].name << " :: " << priors[n].lower << " " << priors[n].upper << endl;
	
	// Run the chains!
	for(int chain = 0; chain < runparams.numchains; chain++){
		runparams.chainID = 1E4 + chain;
		runchain(runparams, data, priors);
	}
	
	cout << "Done" << endl;
		
	
} // END main();

void runchain(struct RPS &runparams, vector<DATA> &data, vector<PARAMP> priors){
	
	// Initialise an MCMC structure
	MarkovChain MCMC;
	// Make sure we start off running the MCMC
	MCMC.run = true;
	// Zero the MCMC step number
	MCMC.step = 0;
	// Zero the acceptance counter
	MCMC.accept_counter = 0;
	// Inherit the number of MCMC steps
	MCMC.numsteps = runparams.numMCMCsteps;
	// Inherit the fraction of MCMC steps to burn.
	MCMC.burninfrac = runparams.burninfrac;
	// How long should the MCMC burn-in time be?
	MCMC.burnintime = floor( MCMC.burninfrac * MCMC.numsteps );

	// Seed the random number generator
	srand (time(NULL));
	
	// The likelihood variables
	double L_current, L_proposed;

	// Get the number of parameters
	int numparams = priors.size();
	
	// Variable holding the values of the parameters
	vector<double> parameters;
	
	// Start off at a random position in parameter space
	for(int param = 0; param < numparams; param++)
		parameters.push_back(priors[param].lower + UnitRand() * (priors[param].upper - priors[param].lower));
	
	// Open up file to dump chain info
	MCMC.chainfile.open( runparams.chaindir + runparams.chainfileprefix + "_" + Int2String(runparams.chainID) + ".dat");
	
	
	while( MCMC.run ){
		
		// Create vector to hold the current values of the parameters
		vector<double> current;
		
		// Dump the values of "parameters" into "current"
		// NOTE: for step = 0, these are the random initial points picked
		for(int n = 0; n < numparams; n ++)
			current.push_back(parameters[n]);
		
		// Current value of the likelihood function
		L_current = ComputeLikelihood(current, data);
		
		// Get some proposed parameters
		vector<double> proposed = GetProposedParameters(priors, current, L_current);
		
		// Get the value of the likelihood with these proposed parameters		
		L_proposed = ComputeLikelihood(proposed, data);
		
		// Check to see if the current step number is within the burn-in time
		if(MCMC.step > MCMC.burnintime){
			
			double LikelihoodRatio = L_proposed / L_current;
			
			// Decide whether to accept the proposed parameters.
			// (1) if the proposed likelihood is larger than current,
			// we accept.
			if( L_proposed >= L_current )
				MCMC.accept = true;
			// However, we also randomly accept according to a 
			// probability:
			else{
				double ran = UnitRand();
				if( ran < LikelihoodRatio ){
					MCMC.accept = true;
				}
				else
					MCMC.accept = false;
			}
			// If we are accepting the proposed values, set the current values
			// to be the proposed ones.
			if(MCMC.accept){
				// Increment acceptance counter
				MCMC.accept_counter++;
				for(int n = 0; n < numparams; n++)
					parameters[n] = proposed[n];
			}
			else{
				for(int n = 0; n < numparams; n++)
					parameters[n] = current[n];
			}

		}
		else{
			// Whilst inside the burn-in time, always update   
			// current values with proposed ones.
			for(int n = 0; n < numparams; n++)
				parameters[n] = proposed[n];
		}
		
		
		// Dump info to file
		for(int n = 0; n < numparams; n++)
			MCMC.chainfile << parameters[n] << " ";
		MCMC.chainfile << " " << L_current << endl;
		
		// Kill MCMC if exceed step number
		if(MCMC.step > MCMC.numsteps) MCMC.run = false;
		MCMC.step ++;				
	}
	
	MCMC.chainfile.close();
		
	cout << "MCMC_accept_counter = " << MCMC.accept_counter << endl;
	
} // END runchain()

double ComputeLikelihood(vector<double> params, vector<DATA> &data){
	
	double test_m = params[0];
	double test_c = params[1];
	double var = 0.0;
	
	for(int n = 0; n < data.size(); n++ ){
		double datax = data[n].x;
		double datay = data[n].y;		
		double testy = test_m * datax + test_c;
		double diff = datay - testy;
		var += diff * diff;
	}
	var /= data.size();
	
	return exp(- 0.5 * var );
	
}// END ComputeLikelihood()



vector<double> GetProposedParameters(vector<PARAMP> priors, vector<double> current, double L_current){
	
	// This is the vector which is to be returned,
	// contains the "successful" proposed parameters (i.e. are inside prior ranges)
	vector<double> proposed;
	
	// Loop over all parameters
	for(int param = 0; param < priors.size(); param++){
		
		// temporary holding variable of the current parameter value
		double thisval = current[param];
		// temporary holding variable of the proposed parameter value
		double prop;
		// Get the upper and lower bounds on the prior for this parameter
		double upper = priors[param].upper;
		double lower = priors[param].lower;
		
		// Get the size of the parameter window allowed by the priors
		double paramwindowsize = upper - lower;
		
		// We will need to keep going untill "GetParamSuccess = true"
		bool GetParamSuccess = false;
		
		// This process may choose parameter values which live outside the prior range,
		// and so must repeat until a parameter is found which is inside
		// prior range.
		while( !GetParamSuccess ){
			
			// Use Box-Muller transform to get N(0,1) -- Normally distributed number.
			prop = thisval + BoxMuller() * min( 1.0 / L_current, paramwindowsize);
			// If the proposed parameter is outside the prior range,
			// we need to do this process again.
			if(prop > upper || prop < lower) 
				GetParamSuccess = false;	
			else
				GetParamSuccess = true;		
		}
		// Dump the sucessful proposed parameter into the vector to be returned.
		proposed.push_back(prop);
	}
	
	// Return the vector containing the proposed parameters back
	return proposed;
	
} // END GetProposedParameters()
