
#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>

#define _USE_MATH_DEFINES

#include "inireader.h"
#include <boost/filesystem.hpp> // Used for creating the output filenames
#include <sstream> // Used for manipulating filenames

void evolve(string *IDS, double *fld, double *runinfo, double *data, double *lagparams, int *impints, int flag);
void GetLagDerivs(double *fld, double *lagparams,int modID, double *lagderivs);
double GetEnergyDensity(double *lagderivs);
double GetPressure(double *lagderivs);
double GetSpeedOfSound(double *lagderivs);
double GetScalarFieldEOM(double *fld,double *lagderivs);
void CheckPathology(double *lagderivs, int flag);
void checkdirexists(string dir);
void writelog(string *IDS, double *runinfo, double *data, int which);
void writeparams(string *IDS, double *runinfo, int nruninfo);
void writescreen(string *IDS, double *runinfo, double *data, int which);

void initialise(int argc, char* argv[],int totID,int ndata,int nruninfo,int nlagparams,int nlagderivs,string *IDS,int *impints,double *fld,double *lagparams,double *runinfo);

using namespace std;
