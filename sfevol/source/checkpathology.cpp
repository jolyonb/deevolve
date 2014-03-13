// checkpathology.cpp

// This checks to see if cs2 < 0 or cs2 > 1.
//  The flag gets altered accordingly.

#include "checkpathology.h"


void CheckPathology(double *lagderivs, int flag){

    double cs2 = GetSpeedOfSound(lagderivs);
    
    if(cs2<0){
        flag=2;
        cout << "Perturbations unstable: cs2 = " << cs2 << endl;
    }
    if(cs2>1){
        flag=3;
        cout << "Superluminal propagation: cs2 = " << cs2 << endl;
    }
    
    
} // END CheckPathology