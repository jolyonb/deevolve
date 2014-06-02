
#include "main.h"


int main(){
	// Start timing!
	boost::timer::cpu_timer myTimer;
	
	// Read in the parameters the user wants to increment
	string UsersInc = "usersinc.ini";
	
	// Output directory for generated .ini files
	string IniGenDir = "ginis/";
	// Prefix for generated ini's
	string GenIniPrefix = "ini_";
	
	// Directory containing the prototypes
	string Pfiledir = "src/prots/";
	// File name of "cosmology" prototype
	string CPFName = "cos.ini";
	// File name of "function" prototype
	string FPFName = "pf.ini";
	// File name of LCDM prototype
	string LCDMfile = "lcdm.ini";
	// File name of "linear-w" prototype
	string LINWfile = "linw.ini";
	// File name of quintessence prototype
	string QUINTfile = "quint.ini";
	// File name of kessence prototype
	string KESSfile = "kess.ini";
	
	// vector to hold users parameter increments
	vector<UIPARS> upars;	
	
	
	// Read in users increment file
	ifstream UI;
	UI.open(UsersInc);

	if(UI){
		while(!UI.eof()){
			UIPARS param;
			UI >> param.name >> param.start >> param.end >> param.inc;
			if(param.name.substr(0,1)!="#") upars.push_back(param);
		}
		UI.close();
	}
	cout << endl;
	cout << "Parameters for incrementation are: " << endl;
	for(int i = 0; i < upars.size(); i++){
		cout << upars[i].name << " " << upars[i].start << " " ;
		cout << upars[i].end << " " << upars[i].inc << endl;
	}
	cout << endl;

	
    IniReader cfile, ffile, mfile;
	ffile.read(Pfiledir+FPFName);
	cfile.read(Pfiledir+CPFName);	
	
	PARAMS params;
	int IDs = 0;
	params.function.ID = IDs; IDs++;
	params.cosmology.ID = IDs;	IDs++;
	params.model.ID = IDs; IDs++;
	
	ReadProts(&params,ffile, params.function.ID);
	ReadProts(&params,cfile, params.cosmology.ID);	
	
	params.model.WhichModel = params.cosmology.model;

	if(params.model.WhichModel == "LambdaCDM")	
		mfile.read(Pfiledir+LCDMfile);
	
	if(params.model.WhichModel == "LinearW")
		mfile.read(Pfiledir+LINWfile);
	
	if(params.model.WhichModel == "Quintessence")
		mfile.read(Pfiledir+QUINTfile);
	
	if(params.model.WhichModel == "Kessence")
		mfile.read(Pfiledir+KESSfile);
	
	
	ReadProts(&params,mfile,params.model.ID);		
	
	
	// Do a quick sanity check
	bool sanity;
	for(int u = 0; u < upars.size(); u++){	
		sanity = false;
		for(int n = 0; n < params.pp.size(); n++){
			if( upars[u].name == params.pp[n].name ){	
				sanity = true;
			}
		}
		if(!sanity) cout << "You picked parameters which dont exist" << endl;
	}	
	
	
	ofstream genparams;
	int fileID = 0;
	string ThisSection = "";
	for(int ngene = 0; ngene < upars.size(); ngene++){
		// Get the incremment on the parameter 
		double pInc = upars[ngene].inc;
		for(double par = upars[ngene].start; par <= upars[ngene].end+pInc; par += pInc, fileID++, genparams.close() ){		
			
			// Open file to dump generated parameters into
			genparams.open(IniGenDir+GenIniPrefix+to_string(10000+fileID)+".gin");
			
			// Dump the new parameter file
			for(int n = 0; n < params.pp.size(); n++){
				// Print section header
				if(params.pp[n].section != ThisSection) {
					genparams << "[" << params.pp[n].section << "]" << endl;
					ThisSection = params.pp[n].section;
				}
				// Print new parameter
				if(	params.pp[n].name == upars[ngene].name )
					genparams << params.pp[n].name << " = " << par << endl;
				else
				// else, print the parameter choices as read in from the prototype ini
			 		genparams << params.pp[n].name << " = " << params.pp[n].val << endl;
			} // END n-loop
		} // END par-loop
	} // END ngene-loop
	
	cout << endl;
	cout << "Finished generating ini files" << endl;
	cout << endl;
    // Stop timing
    myTimer.stop();
    cout << "Generated ini's in " << myTimer.elapsed().wall / 1e6 << " milliseconds." << endl;
	cout << endl;
} // END main()

void ReadProts(struct PARAMS *P, IniReader &init, int ID){	

	if( ID == P->function.ID ){	
		
		PAIRS dp;
		string section = "Function";
		dp.name = "starttime";
		dp.type = "double";
		dp.section = section;
		dp.val = to_string(init.getiniDouble(dp.name, 0.0, section));
		P->pp.push_back(dp);
		
		dp.name = "maxtime";
		dp.type = "double";
		dp.section = section;
		dp.val = to_string(init.getiniDouble(dp.name, 0.0, section));
		P->pp.push_back(dp);
		
		dp.name = "logdir";
		dp.type = "string";
		dp.section = section;
		dp.val =  init.getiniString(dp.name, "logs", section);
		P->pp.push_back(dp);
		
		dp.name = "numberpad";
		dp.type = "int";
		dp.section = section;
		dp.val =  to_string(init.getiniInt(dp.name, 4, section));
		P->pp.push_back(dp);
		
		dp.name = "consistencyclass";
		dp.type = "string";
		dp.section = section;
		dp.val =  init.getiniString(dp.name, "SimpleCheck", section);
		P->pp.push_back(dp);
		
		dp.name = "outputclass";
		dp.type = "string";
		dp.section = section;
		dp.val =  init.getiniString(dp.name, "BasicDump",section);
		P->pp.push_back(dp);
		
		dp.name = "postprocess";
		dp.type = "string";
		dp.section = section;
		dp.val =  init.getiniString(dp.name, "true", section);
		P->pp.push_back(dp);
		
		dp.name = "postname";
		dp.type = "string";
		dp.section = section;
		dp.val =  init.getiniString(dp.name, "d", section);
		P->pp.push_back(dp);
		
		dp.name = "union21";
		dp.type = "string";
		dp.section = section;
		dp.val =  init.getiniString(dp.name, "SCPUnion2.1_mu_vs_z.txt", section);
		P->pp.push_back(dp);
		
	}
	
	if( ID == P->cosmology.ID ){
		
		PAIRS dp;
		string section = "Cosmology";
		
		dp.name = "Tgamma";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "Useh2";
		dp.type = "bool";
		dp.section = section;
		dp.val =  init.getiniString(dp.name, "false", section);
		P->pp.push_back(dp);
		
		dp.name = "Omegamh2";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "Omegabh2";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "Omegakh2";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "Omegam";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "Omegab";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "Omegak";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "Neff";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "phi0";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "phidot0";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "zInit";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "Hubbleh";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "desiredh";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0, section));
		P->pp.push_back(dp);
		
		dp.name = "sigmah";
		dp.type = "double";
		dp.section = section;
		dp.val =  to_string(init.getiniDouble(dp.name, 0,section));
		P->pp.push_back(dp);
		
		dp.name = "model";
		dp.type = "string";
		dp.section = section;
		dp.val =  init.getiniString(dp.name, "MODEL",section);
		P->pp.push_back(dp);
		P->cosmology.model = dp.val;
	
	}
	
	if( ID == P->model.ID){
		
		PAIRS dp;
		string section = P->model.WhichModel;
		
		if(P->model.WhichModel == "LambdaCDM"){
			
			dp.name = "OmegaLambda";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 0.0, section));
			P->pp.push_back(dp);
			
			dp.name = "precise";
			dp.type = "string";
			dp.section = section;
			dp.val =  init.getiniString(dp.name, "true", section);
			P->pp.push_back(dp);
			
		} // END LambdaCDM-ini read
	
		if(P->model.WhichModel == "LinearW"){
			
			dp.name = "wnaught";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 0.0, section));
			P->pp.push_back(dp);
			
			dp.name = "wa";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 0.0, section));
			P->pp.push_back(dp);

			dp.name = "OmegaLambda";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 0.0, section));
			P->pp.push_back(dp);
			
			dp.name = "precise";
			dp.type = "string";
			dp.section = section;
			dp.val =  init.getiniString(dp.name, "true", section);
			P->pp.push_back(dp);

		} // END LinearW-ini read
	
		if(P->model.WhichModel == "Quintessence"){
			
			dp.name = "PotentialType";
			dp.type = "int";
			dp.section = section;
			dp.val = to_string(init.getiniInt(dp.name, 0,section));
			P->pp.push_back(dp);

			dp.name = "mass";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 0.0,section));
			P->pp.push_back(dp);
			
			dp.name = "lambda";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 1.0,section));
			P->pp.push_back(dp);

			dp.name = "alpha";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 1.0,section));
			P->pp.push_back(dp);
			
			dp.name = "beta";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 1.0,section));
			P->pp.push_back(dp);
			
		} // END quintessence-ini read
	
		if(P->model.WhichModel == "Kessence"){
			
			dp.name = "lambda";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 1.0,section));
			P->pp.push_back(dp);

			dp.name = "alpha";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 1.0,section));
			P->pp.push_back(dp);
			
			dp.name = "beta";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 1.0,section));
			P->pp.push_back(dp);
			
			dp.name = "n";
			dp.type = "double";
			dp.section = section;
			dp.val = to_string(init.getiniDouble(dp.name, 1.0,section));
			P->pp.push_back(dp);
			
		} // END Kessence-ini read
		
	}
	
} // END ReadProts()