/*
 * constw.h
 *
 * This implements a constant w dark energy model.
 * Because such a situation can exist in a variety of models, we do not implement the
 * speed of sound or ghost checks.
 *
 */

#ifndef CONSTW_H_
#define CONSTW_H_

#include "model.h"

class ConstW : public Model {

	public:
		// Here are the functions that are overridden by the quintessence class
		int derivatives(const double data[], double derivs[], Parameters &params);
		std::string init(double data[], const double time, Parameters &params, IniReader &init, int &errorstate);

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);

		// Omega_Lambda today
		double OmegaLambda;
		// Equation of state
		double EOSw;

};

#endif /* CONSTW_H_ */
