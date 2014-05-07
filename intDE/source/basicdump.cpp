#include "basicdump.h"

// Function to print data before the run starts
void BasicDump::printinfo(const double data[], IntParams &params) {

	// Require scientific notation, 15 digits
	myLog->precision(15);

	// Print the header
	*myLog << "Omega_m: \t" << params.getparams().OmegaM() << endl;
	*myLog << "Omega_k: \t" << params.getparams().OmegaK() << endl;
	*myLog << "T_gamma: \t" << params.getparams().Tgamma() << endl;
	*myLog << "Omega_r: \t" << params.getparams().OmegaR() << endl;
	*myLog << "h: \t" << params.getparams().h() << endl;
	*myLog << "z_init: \t" << params.getparams().z0() << endl;
	*myLog << "rho_c: \t" << params.getparams().rhoc() << endl << endl;

}

// Function to print a header before output starts (e.g., column headings)
void BasicDump::printheading() {

	// Dumping everything. Make nice headers.
	//*myData << "time, a, redshift, H, Hdot, phi, phidot, phiddot, Omega_m, Omega_r, Omega_k, Omega_Q, w_total, rho_Q/rho_c, P_Q/rho_c, w_Q, Error" << endl;
	*myLog << "Columns in data file are as follows."
		   << endl
		   << "time, a, redshift, H, Hdot, phi, phidot, phiddot, Omega_m, Omega_r, Omega_k, Omega_Q, w_total, rho_Q/rho_c, P_Q/rho_c, w_Q, Error"
		   << endl << endl;

	// Set up the output form
	*myData << scientific;
	myData->precision(15);

}

// Function to print information after each timestep
void BasicDump::printstep(const double data[], const double time, const IntParams &params, const double status[]) {

	// We are just dumping everything here, except for consistency flags
	for (int i = 0; i < 16; i++)
		*myData << status[i] << ", ";
	*myData << status[16] << endl;

}

// Function to print a header before postprocessed output starts (e.g., column headings)
void BasicDump::postprintheading() {

	// Dumping everything. Make nice headers.
	//*myData << "time, a, redshift, H, Hdot, phi, phidot, phiddot, Omega_m, Omega_r, Omega_k, Omega_Q, w_total, rho_Q/rho_c, P_Q/rho_c, w_Q, Error" << endl;
	*myLog << "Columns in post-processed data file are as follows."
		   << endl
		   << "redshift, H, DC, DM, DA, DL, mu"
		   << endl << endl;

	// Set up the output form
	*myPostData << scientific;
	myPostData->precision(15);

}

// Function to print information after each timestep
void BasicDump::postprintstep(const double z, const double H, const double DC, const double DM, const double DA, const double DL, const double mu) {

	// We are just dumping everything here
	*myPostData << z << ", "
			    << H << ", "
			    << DC << ", "
			    << DM << ", "
			    << DA << ", "
			    << DL << ", "
			    << mu << endl;

}

// Constructor
BasicDump::BasicDump(const std::string &filename, const std::string &postname) {

	// Construct filenames
	string logfile = filename + ".log";
	string datfile = filename + ".dat";
	string postdatfile = filename + postname + ".dat";

	// Open the log files
	myLog = new ofstream(logfile.c_str());
	myData = new ofstream(datfile.c_str());
	myPostData = new ofstream(postdatfile.c_str());

}

// Destructor
BasicDump::~BasicDump() {

	// Close the log files...
	myLog->close();
	myData->close();
	myPostData->close();

	// ... and release memory
	delete myLog;
	delete myData;
	delete myPostData;

}

// Did the files open correctly?
bool BasicDump::filesready() {

	if (myLog->is_open() && myData->is_open() && myPostData->is_open())
		return true;
	else
		return false;

}

// Print a "We're done!" message
void BasicDump::printfinish(const double time) {
	*myLog << setprecision(4) << "Evolution complete in " << time << " milliseconds." << endl;
	cout << setprecision(4) << "Evolution complete in " << time << " milliseconds." << endl;
}

// Print a line to the log file
void BasicDump::printlog(const std::string &output) {
	*myLog << output << endl;
}
