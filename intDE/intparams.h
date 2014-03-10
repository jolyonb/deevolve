/*
 * intparams.h
 *
 * This defines the IntParams class, which is essentially a holding class for a Parameters class
 * and a Model class. It is intended that this class is passed into the integration routine.
 *
 */

#ifndef INTPARAMS_H_
#define INTPARAMS_H_

#include <iostream> // for the NULL constant
#include "params.h"
#include "model.h"

class IntParams {

	public:
		// Constructor and destructor
		IntParams(Parameters&, Model&);
		~IntParams();
		// Getters for the two stored classes
		Model &getmodel();
		Parameters &getparams();

	private:
		// Internal storage for the two classes
		Model *myModel;
		Parameters *myParams;

};

#endif /* INTPARAMS_H_ */
