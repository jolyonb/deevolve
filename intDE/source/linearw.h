/*
 * linearw.h
 *
 * This implements a dark energy model with equation of state w = w0 + wa (1 - a).
 * Because such a situation can exist in a variety of models, we do not implement the
 * speed of sound or ghost checks.
 *
 */

#ifndef LINEARW_H_
#define LINEARW_H_

#include "model.h"

class LinearW : public Model {

	public:
		// Here are the functions that are overridden by the LinearW class
		int derivatives(const double data[], double derivs[], Parameters &params);
		std::string init(double data[], const double time, Parameters &params, IniReader &init, int &errorstate);

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);

		// Omega_Lambda today
		double OmegaLambda;
		// Equation of state w(a) = w0 + wa (1 − a)
		double w0;
		double wa;

};

#endif /* LINEARW_H_ */