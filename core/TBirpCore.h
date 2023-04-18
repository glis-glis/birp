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
#include "TData.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TTask.h"
#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorUniform.h"

//--------------------------------------
// TBirpCore
//--------------------------------------

class TBirpCore {
protected:
	std::string _outname;

	// locations, methods, timepoints
	TUniqueContainer<std::string> _locations;
	TData _data;
	TUniqueContainer<size_t> _timepoints;
	TUniqueContainer<std::string> _covariateNames;
	TUniqueContainer<std::string> _speciesNames;

	// read input files
	[[nodiscard]] std::vector<std::string> _getHeader() const;

	[[nodiscard]] bool _checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index,
	                                         const std::string &Title) const;

	void _checkHeader(const coretools::TInputFile &Infile, bool hasCounts, const std::vector<size_t> &countIndices,
	                  const std::vector<size_t> &effortIndices);

	std::vector<size_t> getColumnIndex(const coretools::TInputFile &infile, const std::string &string);
	std::vector<std::string> _getSpeciesOrCovariateNames(const coretools::TInputFile &infile,
	                                                     const std::vector<size_t> &Indices);
	void _readFile(const std::string &Filename, bool hasCounts);

	void _sortTimepoints();

	void _readForInference(bool hasCounts);

	// simulate
	std::vector<size_t> _readN0() const;

	void _fillTimepointsFromCommandLine();

	void _fillLocationsFromCommandLine();

	std::vector<size_t> _simulateSpeciesNames();
	std::vector<size_t> _simulateCovNames();

	// write
	std::vector<std::string> _getMethodspecificHeader(const TMethods &method);
	void _writeFile();

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
