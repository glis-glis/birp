#ifndef TEXAMPLETASK_H_
#define TEXAMPLETASK_H_

#include "TData.h"
#include "TUniqueContainer.h"
#include <string>

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
	TUniqueContainer<std::string> _CIGroupNames;

	// open files
	void _readData();
	static std::vector<std::string> _getAllFilenames(std::string_view Input);

	// read input files
	auto _keepLocations(const std::string &Filename) const;
	void _readFile(const std::string &Filename);
	void _checkIfAllCountsZero();

	void _sortTimepoints();

	std::map<std::string, std::string> _runMCMC(bool Stochastic, const std::map<std::string, std::string> &InitVals);

	// simulate
	void _fillTimepointsFromCommandLine();
	void _fillLocationsFromCommandLine();
	void _fillCIGroupsFromCommandLine();
	std::vector<size_t> _simulateSpeciesNames();
	std::vector<size_t> _simulateCovEffortNames();
	std::vector<size_t> _simulateCovDetectionNames();

	// write
	std::vector<std::string> _getMethodspecificHeader(const TMethods &method);
	void _writeFile(const std::string &Name);
	void _writeTimePoints();
	void _writeCIGroups();
	void _writeInference();
	void _writeMeanParametersState();

	void _clear();

public:
	TBirpCore();

	// main tasks
	void infer();
	void simulate();
};

#endif /* TEXAMPLETASK_H_ */
