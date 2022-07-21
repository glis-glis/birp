/*
 * TBirpCore.cpp
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#include "TBirpCore.h"
#include "TData.h"
#include "TError.h"
#include "TMCMC.h"

using namespace coretools::instances;
using namespace stattools;
using namespace stattools::prior;

//--------------------------------------
// TBirpModel
//--------------------------------------

TBirpModel::TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder)
    // constructor for inference
    : _gamma("gamma", std::make_shared<TUniformFixed<TParameterBase, TypeGamma, 1>>(), {Filename}),
      _effort("effort", std::make_shared<TUniformFixed<TObservationBase, TypeGamma, 1>>(), {Filename}),
      _counts("counts", std::make_shared<TBirpPrior>(&_gamma, &_effort), {}) {

	DAGBuilder.addToDAG(_gamma);
	DAGBuilder.addToDAG(_effort);
	DAGBuilder.addToDAG(_counts);
}

TBirpModel::TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder, size_t Dimensions)
    : TBirpModel(Filename, DAGBuilder) {
	// constructor for simulation
	_counts.getDefinition().setDimensions({Dimensions, Dimensions});
	_counts.getDefinition().setSimulationFile(Filename);
}

//--------------------------------------
// TBirpCore
//--------------------------------------

void TBirpCore::_readCommandLineParameters() {
	using namespace coretools::instances;
	_filename = parameters().getParameterWithDefault<std::string>("out", "birp");
}

void TBirpCore::infer() {
	stattools::TDAGBuilder dagBuilder;

	// read data
	coretools::TMultiDimensionalStorage<TypeCounts, 2> _storageObs;
	TData::read(_storageObs);

	// finalize model
	TBirpModel model(_filename, dagBuilder);

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(_filename, &dagBuilder);
}

void TBirpCore::simulate() {
	stattools::TDAGBuilder dagBuilder;

	// finalize model
	TBirpModel model(_filename, dagBuilder);

	// run MCMC
	stattools::TSimulator simulator;
	simulator.simulate(_filename, &dagBuilder);
}
