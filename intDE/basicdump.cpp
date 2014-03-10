#include "basicdump.h"

// Use the std namespace for the purposes of outputting stuff to screen
using namespace std;

// Function to print data before the run starts
void BasicDump::printinfo(const double data[], IntParams &params) {

	cout << "Beginning run." << endl;
	cout << "Parameters:" << endl;
	cout << "Omega_m: \t" << params.getparams().OmegaM() << endl;
	cout << "Omega_k: \t" << params.getparams().OmegaK() << endl;
	cout << "T_gamma: \t" << params.getparams().Tgamma() << endl;
	cout << "h: \t" << params.getparams().h() << endl;
	cout << "z_init: \t" << params.getparams().z0() << endl;
	cout << endl;

}

// Function to print a header before output starts (e.g., column headings)
void BasicDump::printheading(const double data[], IntParams &params) {

	// We are outputting time, a, phi, phidot, and w_total
	cout << "time \ta \tphi \tphidot \tw_total" << endl;

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

	// We are outputting time, a, phi, phidot, and w_total
	cout << status[0] << "\t" << status[1] << "\t" << status[5] << "\t" << status[6] << "\t" << status[15] << endl;

}
