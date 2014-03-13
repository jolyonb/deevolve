/*
 * quintessence.h
 *
 * This class is an implementation of the Model class; it implements a very basic quintessence model,
 * mostly to make sure that everything is working in the program.
 *
 */

#ifndef QUINTESSENCE_H_
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
		void getstate(const double data[], double, double info[], Parameters &params);
		double speedofsound2(const double data[]);
		bool implementsSOS();
		bool isghost(const double data[]);
		bool implementsghost();
		std::string description();

	private:
		// Here are some internal functions. They're pretty self-explanatory.
		double pressure(const double data[]);
		double energydensity(const double data[]);
		double potential(const double phi);
		double potentialprime(const double phi);
};

#endif /* QUINTESSENCE_H_ */
