
#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES

#include "inireader.h"
#include <boost/filesystem.hpp> // Used for creating the output filenames
#include <boost/timer/timer.hpp>
using namespace std;

#ifndef MAIN_H_
#define MAIN_H_

struct UIPARS{
	string name;
	double start, end, inc;
};

struct FUNCT{
	int ID;	
};

struct COS{
	int ID;
	string model;
};

struct MODEL{
	int ID;
	double OmegaLambda;
	string precise;
	string WhichModel;
};

struct PAIRS{
	string name, val, type, section;
};

struct PARAMS{
	FUNCT function;
	COS cosmology;
	MODEL model;
	vector<PAIRS> pp;
};

#endif

void ReadProts(struct PARAMS *P, IniReader &init, int ID);