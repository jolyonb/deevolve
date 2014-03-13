
#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES

void GetLagDerivs(double *fld, double *lagparams,int modID, double *lagderivs);
double GetEnergyDensity(double *lagderivs);
double GetPressure(double *lagderivs);
double GetSpeedOfSound(double *lagderivs);
double GetScalarFieldEOM(double *fld,double *lagderivs);
void CheckPathology(double *lagderivs, int flag);

using namespace std;