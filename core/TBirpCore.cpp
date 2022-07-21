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
    : _gamma("gamma", std::make_shared<TUniformFixed<TParameterBase, TypeGamma, 1>>(), {Filename}),
      _effort("effort", std::make_shared<TExponentialFixed<TObservationBase, TypeEffort, 2>>(), {}),
      _counts("counts", std::make_shared<TBirpPrior>(&_gamma, &_effort), {}) {

	DAGBuilder.addToDAG(_gamma);
	DAGBuilder.addToDAG(_effort);
	DAGBuilder.addToDAG(_counts);
	DAGBuilder.buildDAG();
}

TBirpModel::TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder,
                       coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
                       coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts)
    // constructor for inference: set storage of counts and efforts
    : TBirpModel(Filename, DAGBuilder) {
	_counts.storage() = Counts;
	_effort.storage() = Efforts;
}

TBirpModel::TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder, size_t NumberOfMethods,
                       size_t NumberOfTimePoints)
    : TBirpModel(Filename, DAGBuilder) {
	// constructor for simulation: set dimensions and file names for simulated data
	_counts.getDefinition().setDimensions({NumberOfMethods, NumberOfTimePoints});
	_counts.getDefinition().setSimulationFile(Filename);

	_effort.getDefinition().setDimensions({NumberOfMethods, NumberOfTimePoints});
	_effort.getDefinition().setSimulationFile(Filename);
}

//--------------------------------------
// TBirpCore
//--------------------------------------

void TBirpCore::_readCommandLineArgs() {
	using namespace coretools::instances;
	_filename = parameters().getParameterWithDefault<std::string>("out", "birp");
}

void TBirpCore::infer() {
	_readCommandLineArgs();

	// read data
	coretools::TMultiDimensionalStorage<TypeCounts, 2> storageCounts;
	coretools::TMultiDimensionalStorage<TypeEffort, 2> storageEfforts;
	TData::read(storageCounts, storageEfforts);

	// create model
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_filename, dagBuilder, storageCounts, storageEfforts);

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(_filename, &dagBuilder);
}

void TBirpCore::simulate() {
	// read command line arguments
	_readCommandLineArgs();
	size_t numTimePoints = parameters().getParameterWithDefault("numTimePoints", 2);
	size_t numMethods    = parameters().getParameterWithDefault("numMethods", 2);

	// create model
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_filename, dagBuilder, numMethods, numTimePoints);

	// run MCMC
	stattools::TSimulator simulator;
	simulator.simulate(_filename, &dagBuilder);
}
