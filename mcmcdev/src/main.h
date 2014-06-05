#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <boost/math/distributions/normal.hpp>
using namespace std;

#define _USE_MATH_DEFINES

#define PI M_PI


double NormalDist(double x, double mu, double sigma){
	
	return exp(-0.5*pow(x-mu,2)*pow(sigma,-2))*pow(2*PI,-0.5)/sigma;
	
} // END NormalDist()

double UniformDist(){
	
	// First, decide what sign the number has:
	
	double sign = rand() / (double)RAND_MAX;
	
	if( sign < 0.5 ) 
		sign = 1.0;
	else
		sign = -1.0;
	
	return sign * rand() / (double)RAND_MAX;
	
} // END UniformInfinite()

double UnitRand(){
	
	return rand() / (double)RAND_MAX;
	
} // END UniformInfinite()

double min(double x1, double x2){
	if(x1 < x2) 
		return x1;
	else
		return x2;
}

double BoxMuller(){
	return sqrt( - 2 * log10( UnitRand() ) ) * cos(2 * PI * UnitRand());
}

double max(double x1, double x2){
	if(x1 > x2) 
		return x1;
	else
		return x2;
}

string Int2String(int Number) {
    return static_cast<ostringstream*>( &(ostringstream() << Number) )->str();
}
