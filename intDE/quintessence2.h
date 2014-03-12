/*
 * quintessence2.h
 *
 * This class is an implementation of the Model class; it implements a very basic quintessence model,
 * mostly to make sure that everything is working in the program.
 *
 * The difference between this quintessence module and quintessence.h is that this module integrates the acceleration equation
 * instead of just using the Friedmann equation to obtain the evolution of a.
 *
 */

#ifndef QUINTESSENCE_H2_
#define QUINTESSENCE_H2_

#include "model.h"
#include "integrate.h"
#include <cmath>

class QuintessenceH : public Model {

	public:
		// Here are the functions that are overridden by the quintessence class
		int derivatives(const double data[], double derivs[], Parameters &params);
		void getstate(const double data[], double, double info[], Parameters &params);
		int init(double data[], double time, Parameters &params);
		double speedofsound2(const double data[], const double derivs[]);
		bool implementsSOS();
		bool isghost(const double data[], const double derivs[]);
		bool implementsghost();

	private:
		// Here are some internal functions. They're pretty self-explanatory.
		double pressure(const double data[]);
		double energydensity(const double data[]);
		double potential(const double phi);
		double potentialprime(const double phi);
};

#endif /* QUINTESSENCE_H2_ */
