#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <string.h>

using namespace std;

string Int2String(int Number) {
    return static_cast<ostringstream*>( &(ostringstream() << Number) )->str();
}

// Structure used to store priors
struct PARAMPRIORS{
	string section;
	string name;
	double lower;
	double upper;
	double sigma;
};
