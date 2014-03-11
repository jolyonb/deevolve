#include "basicdump.h"

// Use the std namespace for the purposes of outputting stuff to screen
using namespace std;

// Function to print data before the run starts
void BasicDump::printinfo(const double data[], IntParams &params) {

	// Require scientific notation, 15 digits
	cout.precision(15);
	cout << scientific;

	// Print the header
	cout << "Beginning run." << endl;
	cout << "Parameters:" << endl;
	cout << "Omega_m: \t" << params.getparams().OmegaM() << endl;
	cout << "Omega_k: \t" << params.getparams().OmegaK() << endl;
	cout << "T_gamma: \t" << params.getparams().Tgamma() << endl;
	cout << "Omega_r: \t" << params.getparams().OmegaR() << endl;
	cout << "h: \t" << params.getparams().h() << endl;
	cout << "z_init: \t" << params.getparams().z0() << endl;
	cout << "rho_c: \t" << params.getparams().rhoc() << endl;
	cout << endl;

}

// Function to print a header before output starts (e.g., column headings)
void BasicDump::printheading(const double data[], IntParams &params) {

	// Dumping everything. Make nice headers.
	cout << "time, a, redshift, H, Hdot, phi, phidot, phiddot, Omega_m, Omega_r, Omega_k, Omega_Q, w_total, rho_Q/rho_c, P_Q/rho_c, w_Q" << endl;

}

// Function to print information after each timestep
void BasicDump::printstep(const double data[], double time, IntParams &params) {

	// An array to hold the status extraction after each step of the integration
	double status[16];

	// Extract the state from the model (not that we presently do anything with it)
	params.getmodel().getstate(data, time, status, params.getparams());

	/* Just as a reminder...
	   * 0 time
	   * 1 a
	   * 2 Redshift
	   * 3 H = \dot{a}/a
	   * 4 \dot{H}
	   * 5 phi
	   * 6 \dot{phi}
	   * 7 \ddot{phi}
	   * 8 Omega_matter (present value)
	   * 9 Omega_radiation (present value)
	   * 10 Omega_k (present value)
	   * 11 Omega_Q (present value)
	   * 12 w_total
	   * 13 rho_Q / rho_c
	   * 14 P_Q / rho_c
	   * 15 w_Q
	 */

	// We are outputting everything here. Intended for dumping to a file and reading in later.
	cout << status[0] << ", ";
	cout << status[1] << ", ";
	cout << status[2] << ", ";
	cout << status[3] << ", ";
	cout << status[4] << ", ";
	cout << status[5] << ", ";
	cout << status[6] << ", ";
	cout << status[7] << ", ";
	cout << status[8] << ", ";
	cout << status[9] << ", ";
	cout << status[10] << ", ";
	cout << status[11] << ", ";
	cout << status[12] << ", ";
	cout << status[13] << ", ";
	cout << status[14] << ", ";
	cout << status[15] << endl;

}
