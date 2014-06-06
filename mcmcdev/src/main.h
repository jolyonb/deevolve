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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;

#define _USE_MATH_DEFINES
const double TwoPi = 2.0 * M_PI;

// Random number generation tools
boost::random::mt19937 RNGtool;
boost::random::uniform_real_distribution<> RNGunit(0.0, 1.0);
boost::random::normal_distribution<> RNGnormal(0.0, 1.0);
// Usage:
// double x = RNGunit(RNGtool) // random real from [0,1)
// double x = RNGnormal(RNGtool) // random real from real distribution with mean 0 and standard deviation 1

double UnitRand(){
	
	//return rand() / (double)RAND_MAX;
    return RNGunit(RNGtool);
	
}

double BoxMuller(){
	
	//return sqrt( - 2 * log( UnitRand() ) ) * cos(TwoPi * UnitRand());
    return RNGnormal(RNGtool);
	
}

string Int2String(int Number) {
	
    return static_cast<ostringstream*>( &(ostringstream() << Number) )->str();
	
}
