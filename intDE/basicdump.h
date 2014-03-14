/*
 * basicdump.h
 *
 * This is an output module. It essentially does the most basic output possible: it dumps everything to screen.
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
		void printfinish(const double time);
		bool filesready();
		void printlog(const std::string&);

		// Constructor
		BasicDump(const std::string &filename = "run");
		// Destructor
		~BasicDump(); // Overrides Output class

	private:
		// File output information
		ofstream *myLog;
		ofstream *myData;

};

#endif /* BASICDUMP_H_ */
