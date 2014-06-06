
#include "chainstats.h"


int main(){
	
	string ChainsDir = "chains/";
	string ChainFileName;
	int numchains = 5;
	int numparams = 2;
	double *Wchain = new double[numparams];
	double *Bchain = new double[];
	int nsamples;
	for(int n = 0; n < numchains; n++){
		ifstream getchain;
		getchain.open(ChainFileName);
		// While this chain is open,
		// Zero the within-chain average
		for(int p=0; p < numparams; p++)
			Wchain[p] = 0.0;
		// Zero the number of samples
		while(!getchain.eof()){
			
			// Read in all the parameters.
			for(int p = 0; p < numparams; p++){
				getchain >> params[p];
				Wchain[p]+=params[p];
			}
			getchain >> likelihood;
			nsamples++;
		}
		for(int p = 0; p < numparams; p++)
			Wchain[p] /= nsamples;
		getchain.close();
		
	}
	
	delete Wchain;
	delete Bchain;
	
	
} // END main()