
// main.h

#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>
#include <vector>
 
#define _USE_MATH_DEFINES

#include "inireader.h"
#include <boost/filesystem.hpp> // Used for creating the output filenames
#include <sstream> // Used for manipulating filenames
void checkdirexists(string dir);
using namespace std;

string GetFileNum(int LogNumber, int padding);

struct UIPARS{
		string name;
		double start, end, inc;
};

struct CHI2PAIRS{
	string name;
	double chi2;
	double chi2max;
	double likelihood;
};
