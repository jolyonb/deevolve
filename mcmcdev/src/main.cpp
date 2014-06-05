
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

struct MarkovChain{
	bool run;
	bool accept;
	int step;
	int numsteps;
	int accept_counter;
	int burnintime;
	double burninfrac;
};


double ComputeLikelihood(vector<double> params, vector<DATA> &data);

vector<double> GetProposedParameters(  vector<PARAMP> priors, vector<double> current, double L_current);

int main(){

	// File name of this chain
	string chaininfofile = "chain.dat";
	// File name of the test data file
	string datafile = "testdata.dat";
	// FIle name of priors file
	string priorsfile = "priors.txt";

	
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
	int numparams = priors.size();
	cout << "The priors are" << endl;
	for(int n = 0; n < numparams; n++)
		cout << priors[n].name << " :: " << priors[n].lower << " " << priors[n].upper << endl;
	
	
	MarkovChain MCMC;
	// Make sure we start off running the MCMC
	MCMC.run = true;
	// 
	MCMC.accept = false;
	// Zero the MCMC step number
	MCMC.step = 0;
	// Zero the acceptance counter
	MCMC.accept_counter = 0;
	// Number of steps to run MCMC for
	MCMC.numsteps = 40000;
	// Fraction of total steps which are to be burnt
	MCMC.burninfrac = 0.1;
	// How long should the MCMC burn-in time be?
	MCMC.burnintime = floor( MCMC.burninfrac * MCMC.numsteps );

	//PFIT fparam;
	srand (time(NULL));
	double L_current;
	double L_proposed;
	
	
	vector<double> parameters;
	for(int param = 0; param < numparams; param++)
		parameters.push_back(priors[param].lower + UnitRand() * (priors[param].upper - priors[param].lower));
	
	// Get a point in the parameter space to start off from
	//m_current = priors[0].lower + UnitRand() * (priors[0].upper - priors[0].lower);
	//c_current = priors[1].lower + UnitRand() * (priors[1].upper - priors[1].lower);
	
	ofstream chaininfo;
	chaininfo.open(chaininfofile);
	
	
	
	while(MCMC.run){
		
		// Create vector to hold the current values of the parameters
		vector<double> current;
		/*
		// Put current values of the parameters into vector
		current.push_back(m_current);
		current.push_back(c_current);
		
		// Current values of the parameters
		fparam.m = m_current;
		fparam.c = c_current;
		*/
		for(int n = 0; n < numparams; n ++)
			current.push_back(parameters[n]);
		
		// Current value of the likelihood function
		L_current = ComputeLikelihood(current, data);
		
		// Get some proposed parameters
		vector<double> proposed = GetProposedParameters(priors, current, L_current);
		/*
		m_proposed = proposed[0];
		c_proposed = proposed[1];
		
		fparam.m = m_proposed;
		fparam.c = c_proposed;
		*/
		// Get the value of the likelihood with these proposed parameters		
		L_proposed = ComputeLikelihood(proposed, data);
		
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
				//m_current = m_proposed;
				//c_current = c_proposed;
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
			//m_current = m_proposed;
			//c_current = c_proposed;	
		}
		
		
		// Dump info to file
		for(int n = 0; n < numparams; n++)
			chaininfo << parameters[n] << " ";
		chaininfo << " " << L_current << endl;
		
		// Kill MCMC if exceed step number
		if(MCMC.step > MCMC.numsteps) MCMC.run = false;
		MCMC.step ++;				
	}
	
	chaininfo.close();
		
		
	cout << "MCMC_accept_counter = " << MCMC.accept_counter << endl;	
	
} // END main();

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
	
	// This is the vector which is to be returned
	vector<double> proposed;
	
	// Loop over all parameters
	for(int param = 0; param < priors.size(); param++){
			
		// temporary holding variable of the proposed parameter value
		double prop;
		// Get the upper and lower bounds on the prior for this parameter
		double upper = priors[param].upper;
		double lower = priors[param].lower;
		double thisval = current[param];
		
		// This process may choose parameter values which live outside the prior range,
		// and so must repeat until a parameter is found which is inside
		// prior range.
		double paramwindowsize = upper - lower;
		bool GetParamSuccess = false;
		while( !GetParamSuccess ){
			
			// Use Box-Muller transform to get x = N(0,1).
			prop = thisval + BoxMuller() * min( 1.0 / L_current, paramwindowsize);
			
			if(prop > upper || prop < lower) 
				GetParamSuccess = false;	
			else
				GetParamSuccess = true;		
		}
		
		proposed.push_back(prop);
	}
	
	return proposed;
	
} // END GetProposedParameters()
