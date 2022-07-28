/*
 * TBirpCore.cpp
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#include "TBirpCore.h"
#include "TMCMC.h"

using namespace coretools::instances;
using namespace stattools;
using namespace stattools::prior;

//--------------------------------------
// TBirpModel
//--------------------------------------

void TBirpModel::_addToDAG(stattools::TDAGBuilder &DAGBuilder) {
	DAGBuilder.addToDAG(&_gamma);
	DAGBuilder.addToDAG(&_effort);
	DAGBuilder.addToDAG(&_counts);
	DAGBuilder.buildDAG();
}

TBirpModel::TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder,
                       coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
                       coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts)
    // constructor for inference: set storage of counts and efforts
    : _gamma("gamma", std::make_shared<TUniformFixed<TParameterBase, TypeGamma, 1>>(), {Filename}),
      _effort("effort", std::make_shared<TUniformFixed<TObservationBase, TypeEffort, 2>>(), Efforts, {}),
      _counts("counts", std::make_shared<TBirpPrior>(&_gamma, &_effort, Filename), Counts, {}) {

	_addToDAG(DAGBuilder);
}

TBirpModel::TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder, size_t NumberOfMethods,
                       size_t NumberOfTimePoints, double LambdaEffort)
    // constructor for simulation: set dimensions and file names for simulated data
    // additionally set exponential prior on efforts
    : _gamma("gamma", std::make_shared<TUniformFixed<TParameterBase, TypeGamma, 1>>(), {Filename}),
      _effort("effort", std::make_shared<TExponentialFixed<TObservationBase, TypeEffort, 2>>(),
              coretools::TMultiDimensionalStorage<TypeEffort, 2>({NumberOfMethods, NumberOfTimePoints}),
              {Filename, coretools::str::toString(LambdaEffort)}),
      _counts("counts", std::make_shared<TBirpPrior>(&_gamma, &_effort, Filename),
              coretools::TMultiDimensionalStorage<TypeCounts, 2>({NumberOfMethods, NumberOfTimePoints}), {Filename}) {

	_addToDAG(DAGBuilder);
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
	double lambdaEffort  = parameters().getParameterWithDefault("lambdaEffort", 1.0);

	// create model
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_filename, dagBuilder, numMethods, numTimePoints, lambdaEffort);

	// run MCMC
	stattools::TSimulator simulator;
	simulator.simulate(_filename, &dagBuilder);
}
