
// main.cpp

#include "main.h"

 
int main(){
    

	// Where are the run*.log files kept?
	string LogsDir = "logs/";
	// Where is the usersinc.ini file?
	string UsersInc = "usersinc.ini";
	// Where to output to
	string dumpdir = "out/";
	// What should the data file with all chi2 be called?
	string dumpfile = "chi2.data";
	checkdirexists(dumpdir);
	
	// Padding on the run-files
	int runfilepadding = 4;
	
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
	// Count the number of logfiles which will need to be read in.
	int NumberOfLogs = 0;
	for(int i = 0; i < upars.size(); i++){
		NumberOfLogs+=(upars[i].end - upars[i].start)/upars[i].inc;
	}
	cout << "Number of logs = " << NumberOfLogs <<endl;

	vector<CHI2PAIRS> WMAP, SN, PLANCK, Hubb, dfGS, SDSS, WiggleZ, Boss9, Boss11;
	double maxW=0, maxS=0, maxP=0, maxH=0, maxd=0, maxSD=0, maxWZ=0, maxB9=0, maxB11=0;
	
	for(int LogNumber = 1; LogNumber <= NumberOfLogs+1; LogNumber++){
		// Open up the logfile
		IniReader logIN;
		
		string ThisLogFile = LogsDir + "run" + GetFileNum(LogNumber, runfilepadding) + ".log";
		logIN.read(ThisLogFile);
		
		/*
		vector<double> items;
		
		// Get all of the parameters we may have incremented,
		// from this logfile.
		for(int n = 0; n < upars.size(); n++){
			double pval = logIN.getiniDouble(upars[n].name, 0.0);
			items.push_back(pval);
		}
		// Dump the parameters which incremented, in order.
		for(int n = 0; n < items.size(); n++){
			OutPut << items[n] << " ";
		}
		*/
		CHI2PAIRS chi2temp;
		
		chi2temp.name = "WMAPchi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxW ) maxW = chi2temp.likelihood;
		WMAP.push_back(chi2temp);
		
		chi2temp.name = "SNchi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxS ) maxS = chi2temp.likelihood;
		SN.push_back(chi2temp);
		
		chi2temp.name = "PLANCKchi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxP ) maxP = chi2temp.likelihood;		
		PLANCK.push_back(chi2temp);
		
		chi2temp.name = "Hubblechi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxH ) maxH = chi2temp.likelihood;
		Hubb.push_back(chi2temp);
		
		chi2temp.name = "6dFGSchi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxd ) maxd = chi2temp.likelihood;
		dfGS.push_back(chi2temp);
		
		chi2temp.name = "SDSSchi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxSD ) maxSD = chi2temp.likelihood;	
		SDSS.push_back(chi2temp);
				
		chi2temp.name = "WiggleZchi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxWZ ) maxWZ = chi2temp.likelihood;		
		WiggleZ.push_back(chi2temp);
		
		chi2temp.name = "BOSSDR9chi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxB9 ) maxB9 = chi2temp.likelihood;	
		Boss9.push_back(chi2temp);
		
		chi2temp.name = "BOSSDR11chi";
		chi2temp.chi2 = logIN.getiniDouble(chi2temp.name, 0);
		chi2temp.likelihood = exp(-0.5*chi2temp.chi2);
		if( chi2temp.likelihood > maxB11 ) maxB11 = chi2temp.likelihood;			
		Boss11.push_back(chi2temp);
		
		
	}
	
	ofstream OutPut;
	OutPut.open(dumpdir+dumpfile);
	
	for(int LogNumber = 1; LogNumber <= NumberOfLogs+1; LogNumber++){
		// Open up the logfile
		IniReader logIN;
		
		string ThisLogFile = LogsDir + "run" + GetFileNum(LogNumber, runfilepadding) + ".log";
		logIN.read(ThisLogFile);
		vector<double> items;
		
		// Get all of the parameters we may have incremented,
		// from this logfile.
		for(int n = 0; n < upars.size(); n++){
			double pval = logIN.getiniDouble(upars[n].name, 0.0);
			items.push_back(pval);
		}
		// Dump the parameters which incremented, in order.
		for(int n = 0; n < items.size(); n++){
			OutPut << items[n] << " ";
		}
		OutPut << WMAP[LogNumber-1].likelihood / maxW << " " ;
		OutPut << SN[LogNumber-1].likelihood / maxS << " " ;
		OutPut << PLANCK[LogNumber-1].likelihood / maxP << " " ;
		OutPut << Hubb[LogNumber-1].likelihood / maxH << " " ;		
		OutPut << dfGS[LogNumber-1].likelihood / maxd << " " ;		
		OutPut << SDSS[LogNumber-1].likelihood / maxSD << " " ;				
		OutPut << WiggleZ[LogNumber-1].likelihood / maxWZ << " " ;						
		OutPut << Boss9[LogNumber-1].likelihood / maxB9 << " " ;								
		OutPut << Boss11[LogNumber-1].likelihood / maxB11 << " " ;										
		OutPut << endl;
	}


	// Close
	OutPut.close();

	cout << "done" << endl;

} // END main

string GetFileNum(int LogNumber, int padding){
	string filenum;
	ostringstream convert;
	convert << LogNumber;
	filenum = convert.str();
	int len = filenum.length();
	for (int i = 0; i < padding - len; i++)
		filenum = "0" + filenum;
	return filenum;
}

void checkdirexists(string dir){
    using namespace boost::filesystem;
    if (!exists(dir + "/")) {
        cout << endl;
        cout << " --> Creating output directory" << endl;
        
        create_directory(dir);
    }
} // END checkdirexists