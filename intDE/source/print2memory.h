/*
 * print2memory.h
 *
 * This defines an output class whereby the log data is saved in memory
 * Note that this does not include the evolution data, just the input parameters and final results
 *
 * This class works by storing information in a property tree. Whenever a value is sent to be printed, it stores it.
 * It also stores some information from the parameters of the evolution.
 *
 * The "printfinal" routine instructs the class to dump the details to file, and then clear the tree.
 * This is useful if one wants to dump multiple runs to one log file.
 *
 * There is an extra routine that allows for properties to be read.
 *
 */

#ifndef PRINT2MEMORY_H_
#define PRINT2MEMORY_H_

#include "output.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using std::string;

class Print2Memory : public Output {
    public:
        // Function to print data before the run starts
        void printinfo(const double data[], Parameters&);

        // Routine to print a key and value to the log file
        // Three overloaded versions
        void printvalue(const std::string&, const std::string&);
        void printvalue(const std::string&, const double);
        void printvalue(const std::string&, const int);

        // Function to print information after everything is finished
        void printfinal(const std::string& = "", const int = 0);


        // Overloaded function to extract saved information
        string getvalue (const string &key, const string &def = "");
        double getvalue (const string &key, const double &def = 0.0);
        double getvalue (const string &key, const int &def = 0);


        // Constructor with optional file name
        Print2Memory(const std::string &filename = "run", const std::string &postname = "d");

        // Virtual destructor
        ~Print2Memory();

        // Check to ensure that output is ready (if not outputting to a file, just return true)
        bool filesready();

    private:
        // The property tree, where everything is stored
        boost::property_tree::ptree data;

        // Output file stream
        std::ofstream *myData;

};

#endif /* PRINT2MEMORY_H_ */
