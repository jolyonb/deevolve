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
#define QUINTESSENCE_H_

#include "model.h"
#include "integrate.h"
#include <cmath>
// stringstreams
#include <iostream>
#include <string>
#include <sstream>

class Quintessence : public Model {

	public:
		// Here are the functions that are overridden by the quintessence class
		int derivatives(const double data[], double derivs[], Parameters &params);
		int init(double data[], double time, Parameters &params, IniReader &init);
		double speedofsound2(const double data[]);
		bool implementsSOS();
		bool isghost(const double data[]);
		bool implementsghost();
		std::string description();

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);
		// These internal functions are specific to quintessence
		double potential(const double phi);
		double potentialprime(const double phi);
};

#endif /* QUINTESSENCE_H2_ */
