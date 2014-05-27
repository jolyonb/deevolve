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
#include "evolve.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>

// Function that finds an appropriate filename (padding is number of characters in the number)
std::string getfilename(const std::string &, const std::string &, const std::string &, const int padding = 4);

#endif /* MAIN_H_ */
