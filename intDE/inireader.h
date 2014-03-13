/*
 * inireader.h
 *
 * This file reads in the initialization file, and provides methods to extract entries from it.
 * Heavily relies upon the BOOST C++ libraries
 *
 */

#ifndef INIREADER_H_
#define INIREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <sstream>

using namespace boost::property_tree;
using namespace std;

class IniReader {
	public:

		void read (const string &filename) {
			// Reads the ini file

			// If no params.ini file exists, use default parameters
			if (!boost::filesystem::exists(filename)) {
				cout << "Warning: " << filename << " does not exist. Using default parameters." << endl;
			} else {
				ini_parser::read_ini(filename, inifile);
			}

		}

		string getiniString (const string &key, const string &def = "", const string &section = "") {
			// Returns a string from the ini file
			// key is the name of the key whose value is desired
			// def is the default to be returned if no value is present
			// section is an optional argument for the section of the ini file (defaults to no section)

			// If we are in a section, construct the tree for that section
			ptree usetree;
			if (section == "")
				usetree = inifile;
			else {
				try {
					usetree = inifile.get_child(section);
				}
				catch (...){
					// Likely got here because that section doesn't exist
					// Just use the default
					return def;
				}
			}

			return usetree.get<string>(key, def);
		}

		double getiniDouble (const string &key, const double &def = 0.0, const string &section = "") {
			// Returns a double from the ini file
			// key is the name of the key whose value is desired
			// def is the default to be returned if no value is present
			// section is an optional argument for the section of the ini file (defaults to no section)

			// If we are in a section, construct the tree for that section
			ptree usetree;
			if (section == "")
				usetree = inifile;
			else {
				try {
					usetree = inifile.get_child(section);
				}
				catch (...){
					// Likely got here because that section doesn't exist
					// Just use the default
					return def;
				}
			}

			return usetree.get<double>(key, def);
		}

		double getiniInt (const string &key, const int &def = 0, const string &section = "") {
			// Returns an integer from the ini file
			// key is the name of the key whose value is desired
			// def is the default to be returned if no value is present
			// section is an optional argument for the section of the ini file (defaults to no section)

			// If we are in a section, construct the tree for that section
			ptree usetree;
			if (section == "")
				usetree = inifile;
			else {
				try {
					usetree = inifile.get_child(section);
				}
				catch (...){
					// Likely got here because that section doesn't exist
					// Just use the default
					return def;
				}
			}

			return usetree.get<int>(key, def);
		}

		double getiniBool (const string &key, const bool &def = false, const string &section = "") {
			// Returns a boolean from the ini file
			// key is the name of the key whose value is desired
			// def is the default to be returned if no value is present
			// section is an optional argument for the section of the ini file (defaults to no section)
			// Note that 1 is true, anything else is false (I think?)

			// If we are in a section, construct the tree for that section
			ptree usetree;
			if (section == "")
				usetree = inifile;
			else {
				try {
					usetree = inifile.get_child(section);
				}
				catch (...){
					// Likely got here because that section doesn't exist
					// Just use the default
					return def;
				}
			}

			return usetree.get<bool>(key, def);
		}

	private:
		// Stores the content of the ini file
		ptree inifile;
};



#endif /* INIREADER_H_ */
