
#include "1Dlike.h"

int main(){
	
	// Which parameter do you want a 1D likelihood plot for?
	string whichparam = "wnaught";
//	string whichparam = "wa";	
	// Where is the priors file?
	string priorsfile = "../priors.txt";
	// How many chains are there?
	int numchains = 5;
	// Where are the chains?
	string chaindir = "../chains/run1/";
	// Whats the prefix on the chain file-name?
	string chainprefix = "chain";
	// Where shall we dump the bins?
	string bindir = "bins/";
	// How many bins in the likelihood?
	int nbins = 25;
	
	cout << "Generating 1D likelihood of " << whichparam << endl;
	
	// Get the priors (to get parameter names & bounds etc)
	ifstream priorsin;
	priorsin.open(priorsfile.c_str());
	vector<PARAMPRIORS> spriors;
	int ID = 0, count = 0;
	if(priorsin){
		while(!priorsin.eof()){
			PARAMPRIORS temp;
			priorsin >> temp.section >> temp.name >> temp.lower >> temp.upper >> temp.sigma;
			spriors.push_back(temp);
			if(temp.name == whichparam) ID = count;
			count ++;
		}
		priorsin.close();
	}
	
	double likelihood;
	int numparams = spriors.size();
	double *parameter = new double[numparams];
	double param_min = spriors[ID].lower;
	double param_max = spriors[ID].upper;

	int *bin = new int[nbins+1];
	double *val = new double[nbins];
	double dparam = (param_max - param_min) / (double)nbins;
	
	// Setup bins
	for(int b = 0; b < nbins+1; b++)
		val[b] = param_min + dparam*b;
	
	// Run over all chains
	int counter;
	for(int n = 0; n < numchains; n++){
		// Zero the bin for this chain
		for(int b = 0; b < nbins; b++) bin[b] = 0;
		counter = 0;
		// Get the file name of this chain
		string thischainfile = chaindir + chainprefix + "_" + Int2String(10000 + n)+".dat";
		ifstream getchain;
		// Open this chain
		getchain.open(thischainfile.c_str());
		if(getchain){
			while(!getchain.eof()){
				// Read in the parameters
				for(int p = 0; p < numparams; p++)
					getchain >> parameter[p];
				
				// Read in the likelihood
				getchain >> likelihood;
				// Increment bin
				for(int b = 0; b < nbins; b++){
					if(parameter[ID] >= val[b] && parameter[ID] < val[b+1]) {
						bin[b]++;
						counter++;
					}
				}
				
			}
			getchain.close();
			// Now dump bin to file			
			ofstream dumpbin;
			string thisbinfile = bindir + whichparam + "_" + Int2String(10000 + n) + ".like";
			dumpbin.open(thisbinfile.c_str());
			for(int b = 0; b < nbins; b++){
				dumpbin << val[b] << "\t" << (double)bin[b]/(double)counter << endl;
			}
			dumpbin.close();	
		}
	}
	
	delete bin;
	delete val;
	
}