/*
 * lambdaCDM.h
 *
 * This is a very basic model - it ignores the scalar field altogether, and implements a cosmological constant with
 * Omega_DE = 0.7 (or as specified).
 *
 */

#ifndef LAMBDACDM_H_
#define LAMBDACDM_H_

#include "model.h"
#include "integrate.h"
#include <cmath>
// stringstreams
#include <iostream>
#include <string>
#include <sstream>

class LambdaCDM : public Model {

	public:
		// Here are the base functions that are overridden by this class
		int derivatives(const double data[], double derivs[], Parameters &params);
		std::string init(double data[], const double time, Parameters &params, IniReader &init, int &errorstate);

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double energydensity(const double data[]);
		double pressure(const double data[], const double hdot);
		// Variable to store Omega_Lambda
		double OmegaLambda;
};

#endif /* LAMBDACDM_H_ */
