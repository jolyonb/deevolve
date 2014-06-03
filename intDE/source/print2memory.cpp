/*
 * print2memory.cpp
 *
 * Contains code supporting the Print2Memory output class
 */

#include "print2memory.h"

// Function to print data before the run starts
void Print2Memory::printinfo(const double data[], Parameters& params){

    // Save the cosmological parameters
    printvalue("OmegaMh2", params.rhoM());
    printvalue("OmegaBh2", params.rhoB());
    printvalue("OmegaKh2", params.rhoK());
    printvalue("OmegaRh2", params.rhoR());
    printvalue("Tgamma", params.Tgamma());
    printvalue("zinit", params.z0());
    printvalue("Neff", params.Neff());

}


// Routine to print a key and value to the log file
// Three overloaded versions
void Print2Memory::printvalue(const std::string& key, const std::string& value){
    data.put(key, value);
}
void Print2Memory::printvalue(const std::string& key, const double value){
    data.put(key, value);
}
void Print2Memory::printvalue(const std::string& key, const int value){
    data.put(key, value);
}

// Function to print information after everything is finished
void Print2Memory::printfinal(const std::string& key, const int ignored){
    // Dump the data of interest to file, and clear the tree

    // Firstly, check to make sure that we actually finished the postprocessing stage
    // If DistanceError or FatalError == 1, then we never performed the postprocessing, so we can bail

    // We're interested in dumping the following information:
    // key, WMAPchi, PLANCKchi, SNchi, Hubblechi, 6dFGSchi, SDSSchi, SDSSRchi, WiggleZchi, BOSSDR9chi, BOSSDR11chi

    *myData << std::scientific << std::setprecision(8) <<
    getvalue(key, 0.0) << "\t" <<
    getvalue("WMAPchi", -1.0) << "\t" <<
    getvalue("PLANCKchi", -1.0) << "\t" <<
    getvalue("SNchi", -1.0) << "\t" <<
    getvalue("Hubblechi", -1.0) << "\t" <<
    getvalue("6dFGSchi", -1.0) << "\t" <<
    getvalue("SDSSchi", -1.0) << "\t" <<
    getvalue("SDSSRchi", -1.0) << "\t" <<
    getvalue("WiggleZchi", -1.0) << "\t" <<
    getvalue("BOSSDR9chi", -1.0) << "\t" <<
    getvalue("BOSSDR11chi", -1.0) << std::endl;

    // Finally, clear the tree
    data.clear();
}


// Overloaded function to extract saved information
string Print2Memory::getvalue (const string &key, const string &def) {
    return data.get<string>(key, def);
}
double Print2Memory::getvalue (const string &key, const double &def) {
    return data.get<double>(key, def);
}
double Print2Memory::getvalue (const string &key, const int &def) {
    return data.get<int>(key, def);
}


// Constructor with optional file name
Print2Memory::Print2Memory(const std::string &filename, const std::string &postname) {

    // Construct filenames
    string datfile = filename + ".dat";

    // Open the log files
    myData = new std::ofstream(datfile.c_str());

}

// Virtual destructor
Print2Memory::~Print2Memory() {

    // Close the log files...
    myData->close();
    // ... and release memory
    delete myData;

}

// Check to ensure that output is ready
bool Print2Memory::filesready() {

    if (myData->is_open())
        return true;
    else
        return false;

}

