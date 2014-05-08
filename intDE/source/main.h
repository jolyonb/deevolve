/*
 * main.h
 *
 * This is just the header file for the main programmatic entry point.
 * It contains references to all the appropriate headers, and defines some functions.
 *
 * This software requires the GSL libraries.
 *
 * This software requires the C++ BOOST libraries (see www.boost.org)
 * These are most easily installed using a package manager (libboost-all-dev on ubuntu)
 *
 */

#ifndef MAIN_H_
#define MAIN_H_

#include "integrate.h"
#include "params.h"
#include "model.h"
#include "quintessence.h"
#include "intparams.h"
#include "output.h"
#include "basicdump.h"
#include "lambdaCDM.h"
#include "consistency.h"
#include "simplecheck.h"
#include "inireader.h"
#include "linearw.h"
#include "kessence.h"
#include "KGB.h"
#include "process.h"

#include <iostream>
#include <cmath>
#include <boost/filesystem.hpp> // Used for creating the output filenames
#include <sstream> // Used for manipulating filenames
#include <algorithm>    // std::reverse
#include <vector>

// Function that controls the evolution
int BeginEvolution(Integrator&, IntParams&, double*, const double, const double, Output&, Consistency&, vector<double>&, vector<double>&);

// Function that the integrator calls to obtain derivatives
int intfunc(double, const double*, double*, void*);

// Function that finds an appropriate filename (padding is number of characters in the number)
std::string getfilename(const std::string &, const std::string &, const std::string &, const int padding = 4);

#endif /* MAIN_H_ */
