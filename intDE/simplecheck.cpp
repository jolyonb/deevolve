#include "simplecheck.h"

using namespace std;

// Checks the state of the system after every timestep
void SimpleCheck::checkstate (double data[], double time, IntParams &params, Output &output, double status[]){
	// We want to check the following conditions:
	// Error in Friedmann equation sufficiently small
	// Equation of state is not phantom
	// Speed of sound is not superluminal
	// Speed of sound is not imaginary
	// Perturbations are not ghost-like

	// Check the error in the Friedmann equation
	if (abs(status[16]) > 1e-9) {
		std::stringstream message;
		message << "Warning: Relative error in Friedmann equation is larger than 1e-9, time t = " << time;
		output.printlog(message.str());
	}

	// Check for phantom equation of state for DE
	if (status[15] < -1) {
		std::stringstream message;
		message << "Warning: Dark energy EOS is phantom, time t = " << time;
		output.printlog(message.str());
	}

	// Check for speed of sound
	if (params.getmodel().implementsSOS()) {
		double speed2 = params.getmodel().speedofsound2(data);
		if (speed2 > 1) {
			std::stringstream message;
			message << "Warning: Speed of sound is superluminal, time t = " << time;
			output.printlog(message.str());
		} else if (speed2 < 0) {
			std::stringstream message;
			message << "Warning: Speed of sound is imaginary, time t = " << time;
			output.printlog(message.str());
		}
	}

	// Check for ghosts
	if (params.getmodel().implementsghost()) {
		if (params.getmodel().isghost(data)) {
			std::stringstream message;
			message << "Warning: Perturbations are ghostlike, time t = " << time;
			output.printlog(message.str());
		}
	}

}

void SimpleCheck::checkfinal (double data[], double time, IntParams &params, Output &output, double status[]){
	// We want to check the following conditions:
	// Equation of state of dark energy is very close to -1
	// Hubble parameter is close to the measured value

	// Check for EOS of DE (warn if the EOS is greater than 0.9)
	if (status[15] > -0.9) {
		std::stringstream message;
		message << "Warning: final equation of state of dark energy is > -0.9";
		output.printlog(message.str());
	}

	// Check for Hubble value
	if (abs(status[3] - 1.0) > 0.1) {
		std::stringstream message;
		message << "Warning: final Hubble parameter is more than 10% away from measured value";
		output.printlog(message.str());
	}

}
