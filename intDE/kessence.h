/*
 * kessence.h
 *
 * This implements a general k-essence class. The user will need to define their own k-essence function.
 * The default implementation is a (\nabla \phi)^4 term.
 *
 */

#ifndef KESSENCE_H_
#define KESSENCE_H_

#include "model.h"

class Kessence : public Model {

	public:
		// Here are the functions that are overridden by the quintessence class
		int derivatives(const double data[], double derivs[], Parameters &params);
		std::string init(double data[], const double time, Parameters &params, IniReader &init, int &errorstate);
		double speedofsound2(const double data[]);
		bool implementsSOS();

		// I'll implement these later
//		bool isghost(const double data[]);
//		bool implementsghost();

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);

		// Just some parameters
		double lambda;
		double alpha;
		double beta;
		double n;

		// Function to calculate the Lagrangian and all its appropriate derivatives:
		// U, Up, Upp, UX, UXX, UXP
		// The results array should be of length 6
		int computelagrangian(const double data[]);
		// Each time the stuff is calculated, store both the data and it, so as not to waste computation time
		double storeddata[4];
		double X;
		double L;
		double LX;
		double LXX;
		double Lp;
		double Lpp;
		double LXp;

};

#endif /* KESSENCE_H_ */
