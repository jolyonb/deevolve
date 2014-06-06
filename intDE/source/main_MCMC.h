
// Header file explicitly for the MCMC code

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

	// Return random number from unit interval
double UnitRand(){
	return RNGunit(RNGtool);
}

// Function to return number from N(0,1): unit normal distribution
double BoxMuller(){	
    return RNGnormal(RNGtool);
}

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

double ComputeLikelihood(IniReader& inifile, string *names, string *sections, double *parameters, int numparams);
void GetProposedParameters(double *priors, double *current, double *proposed, int numparams);
