
// initialise.cpp

// This function reads in parameter files & gets all arrays filled up


#include "main.h"


void initialise(int argc, char* argv[],int totID,int ndata,int nruninfo,int nlagparams,int nlagderivs,string *IDS,int *impints,double *fld,double *lagparams,double *runinfo){


    // Get the parameter file:
    IniReader inifile;
    // If the user specified a param file at runtime
    //  then use that, else, use the default, params.ini
    //   (User specifies via ./EXE my_param.par)
    
    if (argc > 1)
        inifile.read(argv[1]);
    else
        inifile.read("params.ini");
    
    // Read in parameter file (the second entry is a default value if nothing found)
    string  runDIR=inifile.getiniString("runDIR","output");
    string  runIDs=inifile.getiniString("runID","0001");
    
    int     modID=int(inifile.getiniDouble("ModelType",1));
    double  potparam1=inifile.getiniDouble("potparam1",1.0);
    double  potparam2=inifile.getiniDouble("potparam2",1.0);
    double  potparam3=inifile.getiniDouble("potparam3",1.0);
    
    double  z_init=inifile.getiniDouble("z_init",1E4);
    double  phi_init=inifile.getiniDouble("phi_init",2.0);
    double  phidot_init=inifile.getiniDouble("phidot_init",10.0);
    double  OmegaR0=inifile.getiniDouble("OmegaR0",1E-4);
    double  OmegaM0=inifile.getiniDouble("OmegaM0",0.3);
    int     checkpath=int(inifile.getiniDouble("CheckPathology",1));
    
    double  dt=inifile.getiniDouble("TimeStepSize",1E-3);
    int     ttot=int(inifile.getiniDouble("NumTimeSteps",1E9));
    int     verbosity=int(inifile.getiniDouble("Verbosity",1));
    int     filefreq=int(inifile.getiniDouble("FileFreq",1));
    
    string  dsuff=inifile.getiniString("MainDataSuffix","_out.dat");
    string  lsuff=inifile.getiniString("LogFileSuffix","_log.dat");
    string  psuff=inifile.getiniString("ParamsFileSuffix","_params.dat");
    
    
    // Check output directory exists.
    // If not, will be created
    
    checkdirexists(runDIR);
    
    // Start to populate arrays
    
    // Array for strings:

    IDS[0]=runDIR;  // Run Directory
    IDS[1]=runIDs;  // Run ID
    IDS[2]=dsuff;   // Data file suffix
    IDS[3]=lsuff;   // logfile suffix
    IDS[4]=psuff;   // parameter file suffix
    
    // Array for useful integers
    impints[0]=totID;
    impints[1]=ndata;
    impints[2]=nruninfo;
    impints[3]=nlagparams;
    impints[4]=nlagderivs;
    impints[5]=checkpath;
    impints[6]=ttot;
    impints[7]=filefreq;
    impints[8]=verbosity;
    
    // Input field initial conditions
    fld[0] = 1.0 / ( 1.0 + z_init );    // a_init from z_init
    fld[1] = 0.0;                       // H_init (set in main)
    fld[2] = phi_init;
    fld[3] = phidot_init;
    
    // Put parameters used to parameterize the Lagrangian
    lagparams[0] = potparam1;
    lagparams[1] = potparam2;
    lagparams[2] = potparam3;
    
    // Useful runifo
    runinfo[0] = modID;
    runinfo[1] = potparam1;
    runinfo[2] = potparam2;
    runinfo[3] = potparam3;
    runinfo[4] = z_init;
    runinfo[5] = phi_init;
    runinfo[6] = phidot_init;
    runinfo[7] = OmegaR0;
    runinfo[8] = OmegaM0;
    runinfo[9] = dt;
    runinfo[10] = verbosity;
    
} // END initialise


void checkdirexists(string dir){
    
    using namespace boost::filesystem;
    
    cout << endl;
    cout << " --> Creating output directory" << endl;
    
    create_directory(dir);
    
} // END checkdirexists
