// evolve.cpp

// Run over time-steps
// Terminate if a > 1 (sensible...?)
// Also terminate if Omegas chosen badly

#include "evolve.h"


void evolve(string *IDS, double *fld, double *runinfo, double *data, double *lagparams, int *impints, int flag){
    
    
   /*
    impints[0]=totID;
    impints[1]=ndata;
    impints[2]=nruninfo;
    impints[3]=nlagparams;
    impints[4]=nlagderivs;
    impints[5]=checkpath;
    impints[6]=tot;
    */
    
    // Array for RHS of equations of motion
    double  *eom=new double[impints[0]];
    // Array for holding derivatives of the Lagrangian
    double  *lagderivs=new double[impints[4]];
    
    int     totID = impints[0];        // Total number of fields
    int     ttot = impints[6];
    int     ndata = impints[1];
    int     filefreq = impints[7];
    int     checkpath = impints[5];
    int     modID = int(runinfo[0]);
    int     verbosity = impints[8];
    double  OmegaR0 = runinfo[7];
    double  OmegaM0 = runinfo[8];
    double  dt = runinfo[9];
    double  a,H,phi,phidot,a2,H2,tau;
    double  rhoDE,PDE,OmegaM,OmegaR,OmegaDE,wde,wtot,cs2;
    
    
    if(verbosity > 0){
        cout << endl;
        cout << "BEGIN: integrate background field equations" << endl;
        cout << endl;
    }
    
    // Open up output data file
    ofstream dataout;
    string str = IDS[0] + IDS[1] + IDS[2];
	dataout.open(str);
    
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
        GetLagDerivs(fld, lagparams, modID, lagderivs);
        
        // Get rhoDE/rhoc
        rhoDE = GetEnergyDensity( lagderivs );
        // Get PDE/rhoc
        PDE = GetPressure( lagderivs );
        
        // Only on the first time-step, construct H_init from Friedmann eq
        if(t==0){
            fld[1] = sqrt(OmegaM0/a+OmegaR0/a2+a2*rhoDE);
            H = fld[1];
            H2 = pow( H , 2.0 );
            cout << " -> Computed H_init = " << H <<endl;
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
        if(t%filefreq == 0 || t == ttot || t == 0){
            
            // Construct \Omega_X
            
            OmegaR = OmegaR0 / a2 / H2;
            OmegaDE = a2 * rhoDE / H2;
            OmegaM = OmegaM0 / a / H2;
            OmegaM = 1.0 - OmegaR - OmegaDE;
            
            if( OmegaM < 0 || OmegaR < 0 || OmegaDE < 0 ){
                flag = 1;
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
    
    delete  lagderivs;
    delete  eom;
    
} // END evolve

