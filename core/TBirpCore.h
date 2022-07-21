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
#include "TDAGBuilder.h"
#include "TObservationTyped.h"
#include "TParameterTyped.h"
#include "TParameters.h"
#include "TPriorExponential.h"
#include "TPriorUniform.h"
#include "TTask.h"

//------------------------------------------
// TBirpModel
//------------------------------------------

class TBirpModel {
private:
	// all parameters of model (including prior parameters)
	stattools::TParameterTyped<TypeGamma, 1> _gamma;

	// data
	stattools::TObservationTyped<TypeEffort, 2> _effort;
	stattools::TObservationTyped<TypeCounts, 2> _counts;

public:
	TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder);
	// constructor for inference
	TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder,
	           coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
	           coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts);
	// constructor for simulation
	TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder, size_t NumberOfMethods,
	           size_t NumberOfTimePoints);
};

//--------------------------------------
// TBirpCore
//--------------------------------------

class TBirpCore {
private:
	std::string _filename;

	void _readCommandLineArgs();

public:
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
		TBirpCore birp;
		birp.infer();
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
