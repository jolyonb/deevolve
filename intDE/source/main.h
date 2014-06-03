/*
 * main.h
 *
 * This is just the header file for the main programmatic entry point.
 *
 */

#ifndef MAIN_H_
#define MAIN_H_

#include "inireader.h"
#include "params.h"
#include "output.h"
#include "basicdump.h"
#include "print2memory.h"
#include "evolve.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>

// Function that finds an appropriate filename (padding is number of characters in the number)
std::string getfilename(const std::string &, const std::string &, const std::string &, const int padding = 4);

// Routine to perform a single evolution
int doSingleEvolution(IniReader &inifile);

// Routine to perform a sweep over a parameter
int doSweep(IniReader &inifile);

// Structure for storing results from experiments
typedef struct expresults {
    double data[11];
} expresults;

#endif /* MAIN_H_ */
