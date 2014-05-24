/*
 * basicdump.h
 *
 * This is an output module. It dumps basically everything to file, creating a log file that can be machine read as a parameter file.
 *
 */

#ifndef BASICDUMP_H_
#define BASICDUMP_H_

#include "output.h"
#include "intparams.h"
#include <fstream>
#include <boost/timer/timer.hpp>
#include <iomanip>

// Use the std namespace for the purposes of outputting stuff to screen
using namespace std;

class BasicDump : public Output {
	public:
		// Functions that are overridden from the Output class
		void printinfo(const double data[], IntParams&);
		void printheading();
		void printstep(const double data[], const double time, const IntParams&, const double status[]);
		void postprintheading();
		void postprintstep(const double z, const double H, const double DC, const double DM, const double DA, const double DL, const double mu);
		void printfinish(const double time);
		bool filesready();
		void printlog(const std::string&);
		void printvalue(const std::string&, const std::string&);

		// Constructor
		BasicDump(const std::string &filename = "run", const std::string &postname = "d");
		// Destructor
		~BasicDump(); // Overrides Output class

	private:
		// File output information
		ofstream *myLog;
		ofstream *myData;
		ofstream *myPostData;

};

#endif /* BASICDUMP_H_ */
