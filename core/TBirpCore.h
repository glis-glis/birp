#ifndef TEXAMPLETASK_H_
#define TEXAMPLETASK_H_

#include "BirpTypes.h"
#include "TBirpPrior.h"
#include "TData.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TTask.h"
#include "stattools/DAG/TDAGBuilder.h"

//--------------------------------------
// TBirpCore
//--------------------------------------

class TBirpCore {
protected:
	std::string _outname;

	// locations, methods, timepoints
	TUniqueContainer<std::string> _locations;
	TData _data;
	TUniqueContainer<TypeTime> _timepoints;
	TUniqueContainer<std::string> _covariateEffortNames;
	TUniqueContainer<std::string> _covariateDetectionNames;
	TUniqueContainer<std::string> _speciesNames;

	// open files
	static std::vector<std::string> _getAllFilenames();

	// read input files
	auto _keepLocations(const std::string &Filename, bool hasCounts) const;
	void _readFile(const std::string &Filename, bool hasCounts);

	void _sortTimepoints();

	void _readForInference(bool hasCounts);
	void _runMCMC(bool Stochastic);

	// simulate
	void _fillTimepointsFromCommandLine();
	void _fillLocationsFromCommandLine();
	std::vector<size_t> _simulateSpeciesNames();
	std::vector<size_t> _simulateCovEffortNames();
	std::vector<size_t> _simulateCovDetectionNames();

	// write
	std::vector<std::string> _getMethodspecificHeader(const TMethods &method);
	void _writeFile(const std::string &Name);
	void _writeTimePoints();

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
