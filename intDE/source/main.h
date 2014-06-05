/*
 * main.h
 *
 * This is a general header that provides supporting routines to programs implementing the code. As such, this should be considered
 * the one-stop file to include to use all of the available tools in the code.
 *
 * Mostly, it contains includes for the rest of the code, as well as a routine to obtain an appropriate (and available) filename.
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
std::string getfilename(const std::string &dir, const std::string &filebase, const std::string &postbase, const int padding = 4, bool dopostprocess = false) {
    // This routine takes in a directory and an output name
    // It goes and finds the first available filename of the form dir / filebase 00001 etc
    // eg., dir/run00001.log and dir/run00001.dat
    // It checks that both files are free
    // If dopostprocess is true, then also checks for dir/run00001d.dat, where the d is whatever is in postbase
    // Note that even if the output is going to screen, this routine won't make anything bad happen

    using namespace boost::filesystem;

    // Firstly, make sure that the directory exists
    if (!exists(dir + "/")) {
        // Directory doesn't exist. Make it.
        create_directory(dir);
        std::cout << "Creating directory " << dir << "/" << std::endl;
    }

    // Secondly, find a unique filename
    for (int counter = 1; ; counter++) {
        // Construct the file number
        string filenum;
        std::ostringstream convert;
        convert << counter;
        filenum = convert.str();
        // Pad the file number with the appropriate number of zeroes
        int len = filenum.length();
        for (int i = 0; i < padding - len; i++)
            filenum = "0" + filenum;

        // Check for the files
        if (exists(dir + "/" + filebase + filenum + ".log"))
            continue;
        if (exists(dir + "/" + filebase + filenum + ".dat"))
            continue;
        if (dopostprocess && exists(dir + "/" + filebase + filenum + postbase + ".dat"))
            continue;

        // If we got to here, we have a unique filename; return it
        return dir + "/" + filebase + filenum;
    }

    // We really shouldn't get here, but parsers like making sure there's a return
    return "error";

}

#endif /* MAIN_H_ */
