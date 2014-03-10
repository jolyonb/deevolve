#include "intparams.h"

// Constructor. Point the class members to the two objects that were passed in.
IntParams::IntParams(Parameters &params, Model &model){
	myParams = &params;
	myModel = &model;
}

// Remove the references to objects that now no longer need to be pointed to.
IntParams::~IntParams() {
	myParams = NULL;
	myModel = NULL;
}

// Return the model class that was stored
Model &IntParams::getmodel() {
	return *myModel;
}

// Return the parameters class that was stored
Parameters &IntParams::getparams() {
	return *myParams;
}
