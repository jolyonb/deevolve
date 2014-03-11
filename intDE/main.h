/*
 * main.h
 *
 * This is just the header file for the main programmatic entry point. It contains references to all the appropriate headers, and defines some functions.
 *
 */

#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <cmath>
#include "integrate.h"
#include "params.h"
#include "model.h"
#include "quintessence.h"
#include "intparams.h"
#include "output.h"
#include "basicdump.h"
#include "lambdaCDM.h"

// Function that controls the evolution
int BeginEvolution(Integrator&, IntParams&, double*, double, double, Output&);

// Function that the integrator calls to obtain derivatives
int intfunc(double, const double*, double*, void*);

#endif /* MAIN_H_ */
