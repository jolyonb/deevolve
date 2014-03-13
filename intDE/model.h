/*
 * model.h
 *
 * This defines an abstract "Model" class. Any dark energy model can be built off this class. At the very least, it must
 * implement the equations of motion in the derivatives function. Fancier implementations can also compute the speed of sound
 * and check for ghosts.
 *
 * Also available is a state function, which can compute a variety of (useful?) information from the present state.
 *
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "params.h"
#include "inireader.h"
#include <string>

// This defines the Model abstract class.
// Individual models will inherit this class and implement the appropriate functions.
class Model {
	public:
		// The derivatives routine is given the state of the system (a, phi and \dot{\phi} as well as H in some cases)
		//  as well as the parameters of the system, and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi} (and \dot{H}, if necessary))
		virtual int derivatives(const double data[], double derivs[], Parameters &params) = 0;
		// This function returns information on the state of the model
		// A mostly-model independent implementation is provided
		virtual void getstate(const double data[], double time, double info[], Parameters &params);

		// A function to allow the model to initialize itself
		// In particular, this function will have to calculate H from the Friedmann equation as an initial condition
		// if the model is evolving H.
		virtual int init(double data[], double time, Parameters &params, IniReader &init) {return 0;}

		// The speedofsound2 returns the speed of sound squared, given the state of the system
		virtual double speedofsound2(const double data[]) {return 0;}
		// The implementsSOS function returns whether or not a class actually implements the speedofsound2 function
		virtual bool implementsSOS() {return false;}

		// The isghost function is given the state of the system
		// and returns whether or not the theory has become ghostlike
		virtual bool isghost(const double data[]) {return false;}
		// The implementsghost function returns whether or not a class actually implements the isghost function
		virtual bool implementsghost() {return false;}

		// Virtual destructor
		virtual ~Model() {return;}

		// Function to return a description of the model and the parameters its using.
		// Will be outputted to the log file.
		virtual std::string description() {return "";}

	private:
		// Functions to return the energy density (given the data), and the pressure (given data and \dot{H})
		virtual double energydensity(const double data[]) = 0;
		virtual double pressure(const double data[], const double hdot) = 0;

};

#endif /* MODEL_H_ */
