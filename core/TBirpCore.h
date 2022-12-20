/*
 * TBirpCore.h
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#ifndef TEXAMPLETASK_H_
#define TEXAMPLETASK_H_

#include "BirpTypes.h"
#include "TBirpPrior.h"
#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include "coretools/Main/TParameters.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorUniform.h"
#include "coretools/Main/TTask.h"

//------------------------------------------
// TBirpModel
//------------------------------------------

class TBirpModel {
private:
	// all parameters of model (including prior parameters)
	stattools::TParameterTyped<TypeGamma, 1> _gamma;

	void _addToDAG(stattools::TDAGBuilder &DAGBuilder);

public:
	TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder);
};

//--------------------------------------
// TBirpCore
//--------------------------------------

class TBirpCore {
private:
	std::string _outname;

	void _readCommandLineArgs();

public:
	TBirpCore();

	// main tasks
	void infer();
	void simulate();
};

//--------------------------------------
// Tasks
//--------------------------------------

class TTask_infer : public coretools::TTask {
public:
	TTask_infer() { _explanation = "Inferring rates of change in species densities from count data."; };

	void run() override {
		TBirpPrior birp;
		birp.read();
		// birp._writeFile("Test.txt");  // this was to test if file gets saved correctly
		// TBirpCore birp;
	    // birp.infer();
	};
};

class TTask_simulate : public coretools::TTask {
public:
	TTask_simulate() { _explanation = "Simulating count data under Birp model."; };

	void run() override {
		TBirpCore birp;
		birp.simulate();
	};
};

#endif /* TEXAMPLETASK_H_ */
