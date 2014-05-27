/*
 * evolve.h
 *
 * Header file for evolve.cpp, which contains the high level routines to run the evolution.
 *
 */

#ifndef EVOLVE_H_
#define EVOLVE_H_

// Various includes
#include "inireader.h"
#include "params.h"
#include "output.h"
#include "integrate.h"
#include "model.h"
#include "intparams.h"
#include "consistency.h"
#include "simplecheck.h"
#include "process.h"

// Classes for dark energy models
#include "lambdaCDM.h"
#include "linearw.h"
#include "quintessence.h"
#include "kessence.h"
#include "KGB.h"

// Standard C libraries
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>


// Function that performs the evolution
int doEvolution(IniReader&, Parameters&, Output&, vector<double>&, vector<double>&, double&);

// Function that performs post-processing
int PostProcessing(IniReader&, Parameters&, Output&, vector<double>&, vector<double>&);

#endif /* EVOLVE_H_ */