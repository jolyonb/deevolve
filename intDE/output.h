/*
 * output.h
 *
 * This provides an abstract class Output, which output classes inherit.
 *
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "intparams.h"

class Output {
	public:
		// Function to print data before the run starts
		virtual void printinfo(const double data[], IntParams&) = 0;

		// Function to print a header before output starts (e.g., column headings)
		virtual void printheading(const double data[], IntParams&) = 0;

		// Function to print information after each timestep
		virtual void printstep(const double data[], double time, IntParams&) = 0;

		// Virtual destructor
		virtual ~Output() {}

};

#endif /* OUTPUT_H_ */
