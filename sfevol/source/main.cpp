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
    int     nruninfo = 9;     // Number of items constituting information about this run
    int     nlagparams = 3;   // Number of parameters used to parameterize Lagrangian
    int     nlagderivs = 5;   // Number of "derivatives of lagrangian" to be stored
    
    // Array for fields which are updated
    double  *fld=new double[totID];
    // Array for RHS of equations of motion
    double  *eom=new double[totID];
    // Array for data to go into files
    double  *data=new double[ndata];
    // Array for holding info about this run
    double  *runinfo=new double[nruninfo];
    // Array for holding parameters in the Lagrangian
    double  *lagparams=new double[nlagparams];
    // Array for holding derivatives of the Lagrangian
    double  *lagderivs=new double[nlagderivs];
    
    char    outlocedit[50];
    double  H,a,phi,phidot;
    double  OmegaM,OmegaR,OmegaDE ;
    double  tau,a2,H2;
    double  rhoDE,PDE,wde,wtot,cs2;
    int     flag=0;
    string  st2;
    // BEGIN: user input

    
    // Read the params file.
    //  can be specified at runtime via
    //      ./sfevol params1.ini
    //  Default is just params.ini
    
    IniReader inifile;
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
    
    // Construct a_init from z_init
    double  a_init = 1.0 / ( 1.0 + z_init );

    // Construct H_init from Friedmann equation
    // Do this for t=0 only (this happens in main loop)
    double  H_init = 0.0;
    
    // Start to populate arrays
    
    fld[0] = a_init;
    fld[1] = H_init;
    fld[2] = phi_init;
    fld[3] = phidot_init;
    
    lagparams[0] = potparam1;
    lagparams[1] = potparam2;
    lagparams[2] = potparam3;
   
    runinfo[0] = modID;
    runinfo[1] = potparam1;
    runinfo[2] = potparam2;
    runinfo[3] = potparam3;
    runinfo[4] = z_init;
    runinfo[5] = phi_init;
    runinfo[6] =  phidot_init;
    runinfo[7] = OmegaR0;
    runinfo[8] = OmegaM0;
    
    // Print run info to screen
    if(verbosity > 0){
        cout << endl;
        cout << "-----------------------------------------------" << endl;
        cout << "FRW scalar field evolution: conformal time" << endl;
        cout << "J Pearson, Durham, March 2014" << endl;
        cout << "-----------------------------------------------" << endl;
        
        cout << "Run dir: " << runDIR << endl;
        cout << "Run ID: " << runIDs << endl;
        cout << endl;
        cout << "COSMOLOGICAL PARAMETERS: " << endl;
        cout << ":: OmegaM_0 = " << runinfo[8] << endl;
        cout << ":: OmegaR_0 = " << runinfo[7] << endl;
        cout << endl;
        cout << "INITIAL CONDITIONS: " << endl;
        cout << ":: z_init = " << runinfo[4] << endl;
        cout << ":: phi_init = " << runinfo[5] << endl;
        cout << ":: phidot_init = " << runinfo[6] << endl;
        cout << endl;
        cout << "Model info: " << endl;
        cout << "---------------------" << endl;
        cout << "modID = " << runinfo[0];
        
        if(modID==0){
            cout << " :: Quintessence :: phi^2 potential" <<endl;
            cout << "L = X - V(phi), where V(phi) = phi^2 / 2" << endl;
        }
        if(modID==1){
            cout << " :: Quintessence :: exponential potential" <<endl;
            cout << "L = X - V(phi), where V(phi) = e^{- kappa lambda phi}" << endl;
            cout << "lambda = " << runinfo[1] << endl;
        }
        if(modID==2){
            cout << " :: k-essence :: exponential potential" <<endl;
            cout << "L = X^n - V(phi), where V(phi) = e^{-kappa lambda phi}" << endl;
            cout << "lambda = " << runinfo[1] << endl;
            cout << "n = " << runinfo[2] << endl;
        }
        if(modID==3){
            cout << " :: k-essence :: kinetic power-law corrections, exponential potential" <<endl;
            cout << "L = X - V + xi V(X/V)^n, where V(phi) = e^{-kappa lambda phi}" << endl;
            cout << "lambda = " << runinfo[1] << endl;
            cout << "n = " << runinfo[2] << endl;
            cout << "xi = " << runinfo[3] << endl;
        }
    }

    // Open paramsout data file
    ofstream paramout;
    st2=runDIR+runIDs+psuff;
    paramout.open(st2);
    
    for(int ID=0;ID<nruninfo;ID++){paramout << runinfo[ID] << endl;}
    paramout.close();

    // Open log data file
    ofstream logout;
    st2=runDIR+runIDs+lsuff;
    logout.open(st2);
    
    logout << endl;
    logout << "Run location: " << runDIR << endl;
    logout << "Run ID: " << runIDs << endl;
    logout << endl;
    logout << "Model info: " << endl;
    logout << "modID = " << runinfo[0];
    // Write the model specific info
    if(modID==0){
        logout << " :: Quintessence :: phi^2 potential" <<endl;
        logout << "L = X - V(phi), where V(phi) = phi^2 / 2" << endl;
    }
    if(modID==1){
        logout << " :: Quintessence :: exponential potential" <<endl;
        logout << "L = X - V(phi), where V(phi) = e^{- lambda phi}" << endl;
        logout << "lambda = " << runinfo[1] << endl;
    }
    if(modID==2){
        logout << " :: k-essence :: exponential potential" <<endl;
        logout << "L = X^n - V(phi), where V(phi) = e^{- lambda phi}" << endl;
        logout << "lambda = " << runinfo[1] << endl;
        logout << "n = " << runinfo[2] << endl;
    }
    if(modID==3){
        logout << " :: k-essence :: kinetic power-law corrections, exponential potential" <<endl;
        logout << "L = X - V + xi V(X/V)^n, where V(phi) = e^{-kappa lambda phi}" << endl;
        logout << "lambda = " << runinfo[1] << endl;
        logout << "n = " << runinfo[2] << endl;
        logout << "xi = " << runinfo[3] << endl;
    }
    logout << endl;
    logout << " INITIAL CONDITIONS" << endl;
    logout << "---------------------" << endl;
    logout << "z_init = " << runinfo[4] << endl;
    logout << "OmegaR0 = " << runinfo[7] << endl;
    logout << "OmegaM0 = " << runinfo[8] << endl;
    logout << "phi_init = " << runinfo[5] <<endl;
    logout << "phidot_init = " << runinfo[6] << endl;
    logout << endl;
    
    // Dont need the runinfo array any more: delete it
    delete runinfo;
    
    // Open up output data file
    ofstream dataout;  
    st2=runDIR+runIDs+dsuff;
	dataout.open(st2);
    
    // Run over time-steps
    // Terminate if a > 1 (sensible...?)
    // Also terminate if Omegas chosen badly
    if(verbosity>0){
        cout << endl;
        cout << "BEGIN: integrate background field equations" << endl;
        cout << endl;
    }
    
    for(int t=0; t < ttot && a<1 && flag<1 ;t++){

        tau=t*dt;
        
        // Extract fields for ease in the code
        a = fld[0];
        H = fld[1];
        phi = fld[2];
        phidot = fld[3];
        
        // Construct useful auxiliary variables
        a2 = pow( a , 2.0 );
        H2 = pow( H , 2.0 );

        // Model specifics: returns "lagderivs" containing, e.g., {L, LX,...}
        GetLagDerivs(fld,lagparams,modID,lagderivs);
        
        // Get rhoDE/rhoc
        rhoDE = GetEnergyDensity(lagderivs);
        // Get PDE/rhoc
        PDE = GetPressure(lagderivs);
        
        // Only on the first time-step, construct H_init from Friedmann eq
        if(t==0){
            H_init=sqrt(OmegaM0/a+OmegaR0/a2+a2*rhoDE);
            fld[1] = H_init;
            H = fld[1];
            H2 = pow( H , 2.0 );
            cout << " > Computed H_init = " << H_init <<endl;
        }
        
        // Check possible pathology (in sound speed)
        if(checkpath==1){CheckPathology(lagderivs,flag);}
        
        // Construct RHS of equations of motion
        // Equation for adot
        eom[0] = a * H;
        // equation for Hdot
        eom[1] = - 0.5 * H2 - 0.5 * OmegaR0 / a2 - 1.5 * a2 * PDE;
        // equation for phidot
        eom[2] = phidot;
        // equation for phidotdot
        eom[3]= GetScalarFieldEOM(fld,lagderivs);
        
        // Update field values
        for(int ID=0;ID<totID;ID++){fld[ID]=dt*eom[ID]+fld[ID];}
        
        // Write to file every filefreq-timesteps
        if(t%filefreq==0 || t==ttot || t==0){
            
            // Construct \Omega_X
            
            OmegaR = OmegaR0 / a2 / H2;
            OmegaDE = a2 * rhoDE / H2;
            OmegaM = OmegaM0 / a / H2;
            OmegaM = 1.0 - OmegaR - OmegaDE;
            
            if( OmegaM < 0 || OmegaR < 0 || OmegaDE < 0 ){
                flag=1;
                cout << "Pathological Omega_X's" << endl;
                cout << "OmegaR = " << OmegaR << endl;
                cout << "OmegaDE = " << OmegaDE << endl;
                cout << "OmegaM = " << OmegaM << endl;
                cout << endl;
                cout << "timestep number = " << t << endl;
                cout << "a = " << a << endl;
            }
            
            // Construct equations of state
            wtot = ( OmegaR0 / 3.0 + pow( a2, 2.0 ) * PDE ) / ( a * OmegaM0 + OmegaR0 + pow( a2 , 2.0 ) * rhoDE );
            
            wde = PDE / rhoDE;
            
            cs2 = GetSpeedOfSound(lagderivs);
            
            // Dump items into array to write to file...
            data[0] = tau;         // 1. conformal time
            data[1] = fld[0];      // 2. a
            data[2] = fld[1];      // 3. H
            data[4] = fld[2];      // 5. phi
            data[5] = fld[3];      // 6. phidot
            data[6] = OmegaM;      // 7. Omega_M(t)
            data[7] = OmegaR;      // 8. Omega_R(t)
            data[8] = OmegaDE;     // 9. Omega_DE(t)
            data[9] = wtot;        // 10. w_{tot}
            data[10] = wde;        // 11. w_{de}
            data[11] = cs2;        // 12. c_s^2
            
            for(int ID = 0; ID < ndata; ID++){dataout << data[ID] << " " ;}
            dataout << endl;
            
        } // end writetofile


    } // end t-loop
    dataout.close();
    
    /* FOR INFO...
    
    runinfo[0]=modID;
    runinfo[1]=potparam1;
    runinfo[2]=potparam2;
    runinfo[3]=potparam3;
    runinfo[4]=z_init;
    runinfo[5]=phi_init;
    runinfo[6]=phidot_init;
    runinfo[7]=OmegaR0;
    runinfo[8]=OmegaM0;
     
    data[0] = tau;         // 1. conformal time
    data[1] = fld[0];      // 2. a
    data[2] = fld[1];      // 3. H
    data[4] = fld[2];      // 5. phi
    data[5] = fld[3];      // 6. phidot
    data[6] = OmegaM;      // 7. Omega_M(t)
    data[7] = OmegaR;      // 8. Omega_R(t)
    data[8] = OmegaDE;     // 9. Omega_DE(t)
    data[9] = wtot;        // 10. w_{tot}
    data[10] = wde;        // 11. w_{de}
    data[11] = cs2;        // 12. c_s^2

    */
    
    // Write progress to logfile & screen
    
    // If no error: print everything
    if(flag==0){
        logout << endl;
        logout << "Completion without any flagging" << endl;
        logout << endl;
        logout << " FINAL VALUES" << endl;
        logout << "---------------------" << endl;
        logout << "a_final = " << data[1] << endl;
        logout << "z_final = " << 1.0/data[1]-1.0 << endl;
        logout << "OmegaM = " << data[6] << endl;
        logout << "OmegaR = " << data[7] << endl;
        logout << "OmegaDE = " << data[8] << endl;
        logout << "wtot = " << data[9] << endl;
        logout << "wde = " << data[10] << endl;
        logout << "cs2 = " << data[11] << endl;
        logout << endl;
        
        cout << endl;
        cout << "Completion without any flagging" << endl;
        cout << endl;
        cout << "FINAL VALUES" << endl;
        cout << "---------------------" << endl;
        cout << "a_final = " << data[1] << endl;
        cout << "z_final = " << 1.0/data[1]-1.0 << endl;
        cout << "OmegaM = " << data[6] << endl;
        cout << "OmegaR = " << data[7] << endl;
        cout << "OmegaDE = " << data[8] << endl;
        cout << "wtot = " << data[9] << endl;
        cout << "wde = " << data[10] << endl;
        cout << "cs2 = " << data[11] << endl;
        cout << "Done" << endl;
        cout << "---------------" << endl;
        
    }
    // If some error, write info to screen & logfile
    if(flag!=0){
        cout << endl;
        cout << "Ended with problem" << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << "flag = " << flag << endl;
        if(flag==1){cout << "Pathological initial conditions: Omega_tot not unity"<<endl;}
        if(flag==2){cout << "Unstable perturbations: cs2 = " << cs2 << endl; }
        if(flag==3){cout << "Superluminal propagation: cs2 = " << cs2 << endl; }
        cout << endl;
        
        logout << endl;
        logout << "Ended with problem" << endl;
        logout << "!!!!!!!!!!!!!!!!!!!!" << endl;
        logout << "flag = " << flag << endl;
        if(flag==1){logout << "Pathological initial conditions: Omega_tot not unity" << endl;}
        if(flag==2){logout << "Unstable perturbations: cs2 = " << cs2 << endl; }
        if(flag==3){logout << "Superluminal propagation: cs2 = " << cs2 << endl; }
        logout << endl;
    }
    logout.close();
    
    delete  fld;
    delete  eom;
    delete  data;
    delete  lagparams;
    delete  lagderivs;
    
} // end main()

