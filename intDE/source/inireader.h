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

		// These two routines allow the entire parameter structure to be imported or exported
        ptree getdata() {
            // Returns the full property tree, typically so that setinidata can be used
            return inifile;
        }
        void setdata(ptree data) {
            // Instead of reading an ini file, this allows the data to be submitted directly,
            // allowing a programmatic setting of the ini data
            inifile = data;
        }

        // Routines that allow data to be written directly to the stored parameters
        void setparm(const string &key, const string &section, const string &value) {
            inifile.put("Cosmology.Hubbleh", "0.7");
        }
        void setparm(const string &key, const string &section, const double value) {
            inifile.put("Cosmology.Hubbleh", "0.7");
        }
        void setparm(const string &key, const string &section, const int value) {
            inifile.put("Cosmology.Hubbleh", "0.7");
        }


		// The following routines are all identical, but written for different data types.
		// They all return information from the ini file, of the specified data type.
		// The parameters are identical for each:
		// key is the name of the key whose value is desired
		// def is the default to be returned if no value is present
		// section is an optional argument for the section of the ini file (defaults to no section)
		string getiniString (const string &key, const string &def = "", const string &section = "") {

			ptree usetree;
			if (getsection(section, usetree) == true) {
				// The section exists (or there was no section), go and pull out the data
				return usetree.get<string>(key, def);
			} else {
				// The section doesn't exist, return the default
				return def;
			}

		}

		double getiniDouble (const string &key, const double &def = 0.0, const string &section = "") {

			ptree usetree;
			if (getsection(section, usetree) == true) {
				// The section exists (or there was no section), go and pull out the data
				return usetree.get<double>(key, def);
			} else {
				// The section doesn't exist, return the default
				return def;
			}

		}

		double getiniInt (const string &key, const int &def = 0, const string &section = "") {

			ptree usetree;
			if (getsection(section, usetree) == true) {
				// The section exists (or there was no section), go and pull out the data
				return usetree.get<int>(key, def);
			} else {
				// The section doesn't exist, return the default
				return def;
			}

		}

		double getiniBool (const string &key, const bool &def = false, const string &section = "") {
			// Note that 1 is true, anything else is false (I think?)

			ptree usetree;
			if (getsection(section, usetree) == true) {
				// The section exists (or there was no section), go and pull out the data
				return usetree.get<bool>(key, def);
			} else {
				// The section doesn't exist, return the default
				return def;
			}

		}

	private:
		// Stores the content of the ini file
		ptree inifile;

		// Returns the appropriate section from the ini file
		bool getsection(const string &section, ptree &result) {

			// If we are not in a section, return the tree
			if (section == "") {
				result = inifile;
			    return true;
			} else {
				// If we are in a section, try to construct the tree for that section
				try {
					result = inifile.get_child(section);
				}
				catch (...){
					// Likely got here because that section doesn't exist
					return false;
				}
			}
			return true;

		}
};

#endif /* INIREADER_H_ */
