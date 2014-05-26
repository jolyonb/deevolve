/*
 * FR.h
 *
 * This implements a general FR class. The user will need to define their own FR function.
 *
 */

#ifndef FR_H_
#define FR_H_

#include "model.h"

class FR : public Model {

	public:
		// Here are the functions that are overridden by the FR class
		int derivatives(const double data[], double derivs[], Parameters &params);
		std::string init(double data[], const double time, Parameters &params, IniReader &init, int &errorstate);

		// The speed of sound and scalar ghost conditions need to be calculated for FR
		double speedofsound2(const double data[]);
		bool implementsSOS() {return true;}
		bool isghost(const double data[]);
		bool implementsghost() {return true;}
		// The speed of tensor perturbations are unchanged from GR however
		double speedoftensor2(const double data[]) {return 1.0;} 
		bool implementsSOT() {return true;}
		bool istensorghost(const double data[]) {return false;} 
		bool implementstensorghost() {return true;}

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);
		int commonterms(const double data[]);
		// Just some parameters
		double lambda;
		double alpha;
		double beta;
		double n;

		// Function to calculate the Lagrangian and all its appropriate derivatives:
		// The results array should be of length 6
		int computelagrangian(const double data[]);
		// Each time the stuff is calculated, store both the data and it, so as not to waste computation time
		double storeddata[4];
		double a2,a4,a8, hubble;
		double phidot, phidot2, phidot3, phidot4;
		double X;
		double F;
		double Fp;
		double Fpp;
		double Fppp;
};

#endif /* FR_H_ */
