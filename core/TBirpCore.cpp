/*
 * TBirpCore.cpp
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#include "TBirpCore.h"
#include "stattools/MCMC/TMCMC.h"

using namespace coretools::instances;
using namespace stattools;
using namespace stattools::prior;

#include "coretools/Types/strongTypes.h"

//--------------------------------------
// TBirpModel
//--------------------------------------

void TBirpModel::_addToDAG(stattools::TDAGBuilder &DAGBuilder) {
	DAGBuilder.addToDAG(&_gamma);
	DAGBuilder.buildDAG();
}

TBirpModel::TBirpModel(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder)
    : _gamma("gamma", std::make_shared<TUniformFixed<TParameterBase, TypeGamma, 1>>(), {Filename}) {

	_addToDAG(DAGBuilder);
}

//--------------------------------------
// TBirpCore
//--------------------------------------

TBirpCore::TBirpCore(){
	using namespace coretools::instances;
	if(parameters().parameterExists("out")){
		_outname = parameters().getParameter<std::string>("out");
		logfile().list("Will write files with prefix '", _outname, "'. (parameter 'out')");
	} else {
		_outname = "birp";
		logfile().list("Will write files with default prefix '", _outname, "'. (set with 'out')");
	}
}


void TBirpCore::infer() {
	// create model
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_outname, dagBuilder);

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(_outname, &dagBuilder);
}

void TBirpCore::simulate() {
	// create model
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_outname, dagBuilder);

	// run simulator
	stattools::TSimulator simulator;
	simulator.simulate(_outname, &dagBuilder);
}
