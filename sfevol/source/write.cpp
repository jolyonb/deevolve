
// write.cpp

// This contains all the bulky writing routines
// Will write to screen, logfile, & parameters
//  Some are called at the start and end of the code
//      - the argument "which" determines whats printed when.


#include "write.h"

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

void writescreen(string *IDS, double *runinfo, double *data, int which){
    if(runinfo[10]>0){
        if(which==1){
            cout << endl;
            cout << "-----------------------------------------------" << endl;
            cout << "FRW scalar field evolution: conformal time" << endl;
            cout << "J Pearson, Durham, March 2014" << endl;
            cout << "-----------------------------------------------" << endl;
            
            cout << "Run dir: " << IDS[0] << endl;
            cout << "Run ID: " << IDS[1] << endl;
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
            
            if(runinfo[0]==0){
                cout << " :: Quintessence :: phi^2 potential" <<endl;
                cout << "L = X - V(phi), where V(phi) = phi^2 / 2" << endl;
            }
            if(runinfo[0]==1){
                cout << " :: Quintessence :: exponential potential" <<endl;
                cout << "L = X - V(phi), where V(phi) = e^{- kappa lambda phi}" << endl;
                cout << "lambda = " << runinfo[1] << endl;
            }
            if(runinfo[0]==2){
                cout << " :: k-essence :: exponential potential" <<endl;
                cout << "L = X + X^n - V(phi), where V(phi) = e^{-kappa lambda phi}" << endl;
                cout << "lambda = " << runinfo[1] << endl;
                cout << "n = " << runinfo[2] << endl;
            }
            if(runinfo[0]==3){
                cout << " :: k-essence :: kinetic power-law corrections, exponential potential" <<endl;
                cout << "L = X - V + xi V(X/V)^n, where V(phi) = e^{-kappa lambda phi}" << endl;
                cout << "lambda = " << runinfo[1] << endl;
                cout << "n = " << runinfo[2] << endl;
                cout << "xi = " << runinfo[3] << endl;
            }
        } // END which==1
        
        if(which==2){
            
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
            
        } // END which==2
        
        if(which>2){
            cout << endl;
            cout << "Ended with problem" << endl;
            cout << "!!!!!!!!!!!!!!!!!!!!" << endl;
            cout << "flag = " << which-2 << endl;
            if(which-2==1){cout << "Pathological initial conditions: Omega_tot not unity"<<endl;}
            if(which-2==2){cout << "Unstable perturbations: cs2 = " << data[11] << endl; }
            if(which-2==3){cout << "Superluminal propagation: cs2 = " << data[11] << endl; }
            cout << endl;

        }
    }
} // END writescreen

void writeparams(string *IDS, double *runinfo, int nruninfo){
    
    ofstream paramout;
    string str=IDS[0]+IDS[1]+IDS[4];
    paramout.open(str);
    
    for(int ID=0;ID<nruninfo;ID++){
        paramout << runinfo[ID] << endl;
    }
    paramout.close();
    
} // END writeparams

void writelog(string *IDS, double *runinfo, double *data, int which){
    
    ofstream    logout;
    string      st2=IDS[0]+IDS[1]+IDS[3];

    // Write this at the start of a run
    if(which==1){
        

        logout.open(st2);
        logout << endl;
        logout << "Run location: " << IDS[0] << endl;
        logout << "Run ID: " << IDS[1] << endl;
        logout << endl;
        logout << "Model info: " << endl;
        logout << "modID = " << runinfo[0];
        // Write the model specific info
        if(runinfo[0]==0){
            logout << " :: Quintessence :: phi^2 potential" <<endl;
            logout << "L = X - V(phi), where V(phi) = phi^2 / 2" << endl;
        }
        if(runinfo[0]==1){
            logout << " :: Quintessence :: exponential potential" <<endl;
            logout << "L = X - V(phi), where V(phi) = e^{- lambda phi}" << endl;
            logout << "lambda = " << runinfo[1] << endl;
        }
        if(runinfo[0]==2){
            logout << " :: k-essence :: exponential potential" <<endl;
            logout << "L = X^n - V(phi), where V(phi) = e^{- lambda phi}" << endl;
            logout << "lambda = " << runinfo[1] << endl;
            logout << "n = " << runinfo[2] << endl;
        }
        if(runinfo[0]==3){
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
    }


    // Do this at the end of a run
    if(which==2){
        logout.open(st2,std::ofstream::app);
        
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
    }
    
    if(which>2){
        logout.open(st2,std::ofstream::app);
        logout << endl;
        logout << "Ended with problem" << endl;
        logout << "!!!!!!!!!!!!!!!!!!!!" << endl;
        logout << "flag = " << which-2 << endl;
        if(which-2==1){logout << "Pathological initial conditions: Omega_tot not unity" << endl;}
        if(which-2==2){logout << "Unstable perturbations: cs2 = " << data[11] << endl; }
        if(which-2==3){logout << "Superluminal propagation: cs2 = " << data[11] << endl; }
        logout << endl;
        
    }

    logout.close();
    
} // END writelog


