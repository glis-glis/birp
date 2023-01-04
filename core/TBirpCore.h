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
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TTask.h"
#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorUniform.h"

//------------------------------------------
// TBirpModel
//------------------------------------------

class TBirpModelPoisson {
private:
	// all parameters of model (including prior parameters)
	stattools::TParameterTyped<TypeGamma, 1> _gamma;
	stattools::TParameterTyped<TypeTime, 1> _timesOfChange;

	std::unique_ptr<stattools::TObservationTyped<TypeCounts, 2>> _data;

public:
	TBirpModelPoisson(std::vector<TSpecies> *Species, TUniqueContainer<std::string> *Locations,
	                  TUniqueContainer<std::string> *Methods, TUniqueContainer<std::size_t> *TimePoints,
	                  const std::string &Filename, stattools::TDAGBuilder &DAGBuilder);
};

//--------------------------------------
// TBirpCore
//--------------------------------------

class TBirpCore {
private:
	std::string _outname;

	// data
	std::vector<TSpecies> _species;

	// locations, methods, timepoints
	TUniqueContainer<std::string> _locations;
	TUniqueContainer<std::string> _methods;
	TUniqueContainer<size_t> _timepoints;

	// read input files
	std::vector<std::string> _getHeader() const;
	void _checkHeader(const coretools::TInputFile &Infile);
	bool _checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index, const std::string &Title) const;
	void _readFile(const std::string &Filename, bool hasCounts);
	void _sortTimepoints();
	void _readForInference();

	// simulate
	std::vector<size_t> _readN0() const;
	void _fillTimepointsFromCommandLine();
	void _fillLocationsFromCommandLine();
	void _fillMethodsFromCommandLine();

	// write
	void _writeFile(const std::string &Filename);

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
