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

class LambdaCDM : public Model {

	public:
		// Here are the functions that are overridden by the quintessence class
		int derivatives(const double data[], double derivs[], Parameters &params);
		void getstate(const double data[], double, double info[], Parameters &params);

	private:
		// Here are some internal functions. They're pretty self-explanatory.
		double pressure(const double data[]);
		double energydensity(const double data[]);
};

#endif /* LAMBDACDM_H_ */
