// model.cpp

// This file contains everything to do with the scalar field model.

// The functions here are

//  (1) GetLagDerivs()
//  (2) GetEnergyDensity()
//  (3) GetPressure()
//  (4) GetSpeedOfSound()
//  (5) GetScalarFieldEOM()

// (1) GetLagDerivs()
// INPUT:
//      fld:        field array
//      lagparams:  quantities used to parameterise the Lagrangian
//      modID:      id number of the model
// OUTPUT:
//      lagderivs:  array containing all useful derivatives of the Lagrangian

// (2) GetEnergyDensity()
// INPUT:   lagderivs
// OUTPUT:  double, rho_DE/rho_c = ( 2 X LX - L ) / 3

// (3) GetPressure()
// INPUT:   lagderivs
// OUTPUT:  double, P_DE/rho_c = L/3

// (4) GetSpeedOfSound()
// INPUT:   lagderivs
// OUTPUT:  double, cs2 = LX / (LX + 2 X LXX )

// (5) GetScalarFieldEOM()
// INPUT:   fld, lagderivs
// OUTPUT:  double, phidotdot

#include "main.h"

void GetLagDerivs(double *fld, double *lagparams, int modID, double *lagderivs){
    
    // Returns derivatives of the Lagrangian
    //  This is the only place that needs to be modified to
    //   modify for your parameterization of the Lagrangian.
    
    // At the moment, k-essence only
    
    double  a = fld[0];
    double  phi = fld[2];
    double  phidot = fld[3];
    
    double  X = pow( phidot , 2.0 ) / a / a / 2.0;
    
    double  potparam1=lagparams[0];
    double  potparam2=lagparams[1];
    double  potparam3=lagparams[2];
    
    double  L,LX,LXX,Lp,LpX,V,dV;
    
    // phi^2 quintessence
    if( modID == 0 ){ 
        
        V = pow( phi , 2.0 ) / 2.0;
        dV = phi;
        L = X - V;      // L
        LX = 1.0;       // dL/dX
        LXX = 0.0;      // d2L/dX2
        Lp = - dV;      // dL/dphi
        LpX = 0.0;      // d2L/dXdphi
        
    } // end modID==1

    
    // exponential quintessence
    if( modID == 1 ){
        
        V = exp( - phi * potparam1 );
        dV = - potparam1 * V;
        L = X - V;      // L
        LX = 1.0;       // dL/dX
        LXX = 0.0;      // d2L/dX2
        Lp = - dV;      // dL/dphi
        LpX = 0.0;      // d2L/dXdphi
        
    } // end modID==1
    
    // power-law k-essence with exp potential
    if( modID == 2 ){
        
        V = exp( - phi * potparam1 );
        dV = - potparam1 * V;
        
        L = pow( X , potparam2 ) - V;
        LX = potparam2 * pow( X , potparam2 - 1.0 );
        LXX = potparam2 * ( potparam2 - 1.0 ) * pow( X , potparam2 - 2.0 );
        Lp = - dV;
        LpX = 0.0;
        
    } // end modID==2
    
    // see Tamanini paper
    if( modID == 3 ){
        
        // L = X - V + \xi V (X/V)^n
        
        // V = e^{-kappa phi lambda}
        // xi = potparam3
        // n = potparam2
        // lambda = potparam3
        
        V=exp(-phi*potparam1);
        dV=-potparam1*V;
        
        L=X-V+potparam3*V*pow(X/V,potparam2);      // L
        LX=1.0+potparam3*potparam2*pow(X,potparam2-1.0)*pow(V,1.0-potparam2);       // dL/dX
        LXX=potparam3*potparam2*(potparam2-1.0)*pow(X,potparam2-2.0)*pow(V,1.0-potparam2);      // d2L/dX2
        Lp=-dV+potparam3*(1.0-potparam2)*pow(X,potparam2)*pow(V,-potparam2);      // dL/dphi
        LpX=potparam3*potparam2*(1.0-potparam2)*pow(X,potparam2-1.0)*pow(V,-potparam2);      // d2L/dXdphi
        
    } // end modID==3
    
    
    lagderivs[0]=X;
    lagderivs[1]=L;
    lagderivs[2]=LX;
    lagderivs[3]=LXX;
    lagderivs[4]=Lp;
    lagderivs[5]=LpX;
    
} // END GetLagDerivs

double GetEnergyDensity(double *lagderivs){
    
    // Returns rho_DE/rho_c
    
    double  X=lagderivs[0];
    double  L=lagderivs[1];
    double  LX=lagderivs[2];
    
    return ( 2.0 * X * LX - L ) / 3.0;
    
}// END GetEnenergyDensity

double GetPressure(double *lagderivs){
    
    // Returns P_DE/rho_c
    
    double  L=lagderivs[1];
    
    return L / 3.0;
    
}// END GetPressure

double GetSpeedOfSound(double *lagderivs){
    
    // Returns scalar field speed of sound
    
    double  X = lagderivs[0];
    double  LX = lagderivs[2];
    double  LXX = lagderivs[3];
    
    return LX / ( LX + 2.0 * X * LXX );
    
} // END GetSpeedOfSound

double GetScalarFieldEOM(double *fld,double *lagderivs){
    
    // Returns scalar field EoM: phidotdot
    
    double  a=fld[0];
    double  H=fld[1];
    double  phidot=fld[3];
    
    double  a2=pow(a,2.0);
    double  LX=lagderivs[2];
    double  Lp=lagderivs[4];
    double  LpX=lagderivs[5];
    double  cs2=GetSpeedOfSound(lagderivs);
    
    return ( 1.0 - 3.0 * cs2 ) * phidot * H + cs2 / LX * ( a2 * Lp - pow( phidot , 2.0 ) * LpX );
    
    
} // END GetScalarFieldEOM

