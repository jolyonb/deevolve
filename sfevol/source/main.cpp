// main.cpp

////////////////////////////////////////////////////////////////////
//
//      Scalar fields in conformal time: evolver
//          J. Pearson, Durham, March 2014
//
//      Main entry into code
//
//      This file contains mainly declarations & output statements
//
////////////////////////////////////////////////////////////////////

#include "main.h"

int main(int argc, char* argv[]){

    int     totID = 4;        // Total number of fields
    int     ndata = 12;       // number of items dumped into file
    int     nruninfo = 11;     // Number of items constituting information about this run
    int     nlagparams = 3;   // Number of parameters used to parameterize Lagrangian
    int     nlagderivs = 5;   // Number of "derivatives of lagrangian" to be stored
    
    // Array for fields which are updated
    double  *fld=new double[totID];

    // Array for data to go into files
    double  *data=new double[ndata];
    // Array for holding info about this run
    double  *runinfo=new double[nruninfo];
    // Array for holding parameters in the Lagrangian
    double  *lagparams=new double[nlagparams];
    // Array for important integers
    int     *impints=new int[9];
    // Array of strings: will hold runDIR & runID
    string *IDS=new string[5];
    
    // Warning flag: this should remain=0 throughout
    //  evolution, unless something goes wrong.
    //      (output of logfile will tell you what happend)
    int     flag = 0;
  
    // Go read params file, and fill up arrays
    initialise(argc,argv,totID,ndata,nruninfo,nlagparams,nlagderivs,IDS,impints,fld,lagparams,runinfo);

    // Write initial to screen
    writescreen(IDS,runinfo,data,1);

    // Write parameters to file
    writeparams(IDS,runinfo,nruninfo);
    // Write logfile
    writelog(IDS,runinfo,data,1);
    // Evolve
    evolve(IDS,fld,runinfo,data,lagparams,impints,flag);
    
    // Write final info to screen
    writescreen(IDS,runinfo,data,2+flag);
    // Write final info to logfile
    writelog(IDS,runinfo,data,2+flag);
        
    delete  fld;
    delete  data;
    delete  lagparams;
    delete  impints;
    delete  runinfo;

    
} // end main()

