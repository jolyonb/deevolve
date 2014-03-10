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

class BasicDump : public Output {
	public:
		// Functions that are overridden from the Output class
		void printinfo(const double data[], IntParams&);
		void printheading(const double data[], IntParams&);
		void printstep(const double data[], double time, IntParams&);
};

#endif /* BASICDUMP_H_ */
