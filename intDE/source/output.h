/*
 * output.h
 *
 * This provides an abstract class Output, which output classes inherit.
 *
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iostream>
#include "intparams.h"

class Output {
	public:
		// Function to print data before the run starts
		virtual void printinfo(const double data[], IntParams&) = 0;

		// Function to print a header before output starts (e.g., column headings)
		virtual void printheading() {}

		// Function to print information after each timestep
		virtual void printstep(const double data[], const double time, const IntParams&, const double status[]) = 0;

		// Function to print header for postprocessing data before output starts (e.g., column headings)
		virtual void postprintheading() {}

		// Function to print information for postprocessing data after each time step
		virtual void postprintstep(const double z, const double H, const double DC, const double DM, const double DA, const double DL, const double mu) = 0;

		// Function to print information after run is complete (time is given in milliseconds)
		virtual void printfinish(const double time) {}

		// Routine to print a line to the log file
		virtual void printlog(const std::string&) = 0;

		// Routine to print a value and line to the log file
		virtual void printvalue(const std::string&, const std::string&) = 0;

		// Constructor with optional file name
		Output(const std::string &filename = "run", const std::string &postname = "d") {}

		// Virtual destructor
		virtual ~Output() {}

		// Check to ensure that output is ready (if not outputting to a file, just return true)
		virtual bool filesready() = 0;

};

#endif /* OUTPUT_H_ */
