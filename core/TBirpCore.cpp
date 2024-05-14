/*
 * TBirpCore.cpp
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#include "TBirpCore.h"
#include "TData.h"
#include "coretools/Types/strongTypes.h"
#include "coretools/devtools.h"
#include "stattools/MCMC/TMCMC.h"
#include "stattools/Priors/TPriorNormal.h"
#include <set>

using namespace coretools::instances;
using namespace stattools;
using namespace stattools::prior;

//--------------------------------------
// TBirpCore
//--------------------------------------

TBirpCore::TBirpCore() {
	using namespace coretools::instances;
	if (parameters().exists("out")) {
		_outname = parameters().get<std::string>("out");
		logfile().list("Will write files with prefix '", _outname, "'. (argument 'out')");
	} else {
		_outname = "birp";
		logfile().list("Will write files with default prefix '", _outname, "' (use argument 'out' to change).");
	}
}

void TBirpCore::_readForInference(bool hasCounts) {
	// get all input filenames (one per method)
	std::vector<std::string> filenames = _getAllFilenames();

	// read files (one per method)
	for (const auto &name : filenames) { _readFile(name, hasCounts); }

	// sort and finalize
	_sortTimepoints();
	_data.fillMethLocIndices(_speciesNames.size(), _locations);

	// write filtered counts
	_writeFile("filtered");
	_writeTimePoints();

	logfile().endIndent();
}

//-----------------------------------
// Get all filenames

std::string getDelim(std::string_view Filename) {
	std::string delim       = "\t";
	std::string_view ending = coretools::str::readAfterLast(Filename, '.');
	if (ending == "csv") { delim = ","; }
	return delim;
}

std::vector<std::string> readFilenamesFromFile(std::string_view Filename) {
	// open file. Can be either:
	// 1) the file containing the counts for a single method
	// 2) a file containing all filenames for multiple methods
	coretools::TInputMaybeRcppFile file(Filename, coretools::FileType::Header, getDelim(Filename));

	if (file.header().size() != 1) { // case 1: directly counts file
		return {(std::string)Filename};
	}
	// case 2: all filenames
	// open again without header
	coretools::TInputMaybeRcppFile file2(Filename, coretools::FileType::NoHeader, getDelim(Filename));
	std::vector<std::string> filenames;
	for (; !file2.empty(); file2.popFront()) { filenames.emplace_back(file2.get(0)); }
	return filenames;
}

std::vector<std::string> TBirpCore::_getAllFilenames() {
	// read --data -> this can be:
	// 1) the file containing the counts for a single method
	// 2) a file containing all filenames for multiple methods
	// 3) a comma-separated list of filenames for multiple methods
	std::string input = parameters().get("data");
	logfile().startIndent("Will read input data from ", input, " (argument 'data').");

	std::vector<std::string> filenames;
	if (std::filesystem::exists(input)) { // is a file
		filenames = readFilenamesFromFile(input);
	} else {
		coretools::str::fillContainerFromString(input, filenames, ',');
	}

	if (filenames.empty()) { UERROR("Provided filenames are empty (argument 'data')."); }

	return filenames;
}

//-----------------------------------
// Read one method's file

std::vector<std::string> getHeader() { return {"location", "timepoint"}; }

void checkHeader(const coretools::TInputMaybeRcppFile &Infile, bool hasCounts, const std::vector<size_t> &countIndices,
                 const std::vector<size_t> &CovEffortIndices) {
	std::vector<std::string> headerVec = getHeader();

	// make sure location and timepoint are included in header
	for (const auto &colName : getHeader()) {
		if (std::find(Infile.header().begin(), Infile.header().end(), colName) == Infile.header().end()) {
			UERROR("Error in header of file ", Infile.name(), ": Mandatory column '", colName, "' is missing.");
		}
	}

	// make sure there is at least one column starting with "covEffort" and counts
	if (CovEffortIndices.empty()) {
		UERROR("Error in header of file ", Infile.name(), ": Need at least one column starting with 'covEffort_'.");
	}
	if (hasCounts && countIndices.empty()) {
		UERROR("Error in header of file ", Infile.name(), ": Need at least one column starting with 'counts_'.");
	}
}

std::vector<size_t> getColumnIndex(const coretools::TInputMaybeRcppFile &infile, const std::string &string) {
	std::vector<size_t> ix;
	for (size_t i = 0; i < infile.numCols(); i++) {
		if (coretools::str::stringStartsWith(infile.header()[i], string)) { ix.push_back(i); }
	}
	return ix;
}

std::vector<std::string> getSpeciesOrCovariateNames(const coretools::TInputMaybeRcppFile &infile,
                                                    const std::vector<size_t> &Indices) {
	std::vector<std::string> names(Indices.size());
	for (size_t s = 0; s < Indices.size(); s++) {
		names[s] = coretools::str::readAfter(infile.header()[Indices[s]], '_');
	}
	return names;
}

std::vector<size_t> fillCovNames(coretools::TInputMaybeRcppFile &infile, const std::vector<size_t> &Indices,
                                 TUniqueContainer<std::string> &UniqueContainer) {
	std::vector<std::string> covNames = getSpeciesOrCovariateNames(infile, Indices);
	std::vector<size_t> cov_ids(covNames.size());
	for (size_t c = 0; c < covNames.size(); ++c) { cov_ids[c] = UniqueContainer.add(covNames[c]); }
	return cov_ids;
}

std::string getMethodName(const std::string &Filename) {
	std::string_view methodNameFull = Filename;
	if (coretools::str::stringContains(Filename, '/')) {
		methodNameFull = coretools::str::readAfterLast(Filename, '/');
	}
	std::string methodName = (std::string)coretools::str::readBeforeLast(methodNameFull, '.');
	coretools::str::eraseAllOccurences(methodName, "_simulatedCounts");
	return methodName;
}

auto getIndicesMethodFile(const coretools::TInputMaybeRcppFile &Infile, bool hasCounts) {
	// get position of columns that contain counts: might be an empty vector in case the file has no counts
	std::vector<size_t> countIndices = getColumnIndex(Infile, "counts");

	// get position of columns that contain covariates (efforts)
	std::vector<size_t> covEffortIndices    = getColumnIndex(Infile, "covEffort");
	std::vector<size_t> covDetectionIndices = getColumnIndex(Infile, "covDetection");

	checkHeader(Infile, hasCounts, countIndices, covEffortIndices);

	auto methodName = getMethodName(Infile.name());
	return std::make_tuple(countIndices, covEffortIndices, covDetectionIndices, methodName);
}

bool allEffortsAreZero(coretools::TInputMaybeRcppFile &InputFile, const std::vector<size_t> &CovEffortIndices) {
	std::set<TypeCovariateEffort> covEffort;
	for (auto c : CovEffortIndices) { covEffort.insert(InputFile.get<TypeCovariateEffort>(c)); }
	return covEffort.size() == 1 && *covEffort.begin() == 0.0; // returns true if all covariates of effort are zero
}

auto TBirpCore::_keepLocations(const std::string &Filename, bool hasCounts) const {
	// open file and get relevant column indices
	coretools::TInputMaybeRcppFile file(Filename, coretools::FileType::Header, getDelim(Filename));
	auto [countIndices, covEffortIndices, covDetectionIndices, methodName] = getIndicesMethodFile(file, hasCounts);

	// map: per location, how many timepoints with non-zero effort are there?
	std::map<std::string, size_t> numTimePointsPerLocation;
	for (; !file.empty(); file.popFront()) {
		// read location and timepoint
		std::string location  = file.get<std::string>("location");
		std::string timepoint = file.get<std::string>("timepoint");

		// read effort(s)
		if (allEffortsAreZero(file, covEffortIndices)) {
			// all effort covariates are zero
			logfile().warning("Method ", methodName, ": Ignoring location ", location, " at timepoint ", timepoint,
			                  " because of zero effort.");
			continue;
		}

		auto it = numTimePointsPerLocation.find(location);
		if (it == numTimePointsPerLocation.end()) { // new location
			numTimePointsPerLocation[location] = 1;
		} else {
			it->second++; // location already exists
		}
	}

	// parsed entire file: now check for each location if there are enough timepoints
	std::map<std::string, bool> keepLocations;
	size_t numLocationsToKeep = 0;
	for (const auto &loc : numTimePointsPerLocation) {
		if (loc.second == 1) {
			logfile().warning("Method ", methodName, ": Ignoring location ", loc.first,
			                  " because it only has a single timepoint with non-zero effort.");
			keepLocations[loc.first] = false;
		} else {
			keepLocations[loc.first] = true;
			++numLocationsToKeep;
		}
	}
	bool keepMethod = true;
	if (numLocationsToKeep == 0) {
		logfile().warning("Ignoring method ", methodName,
		                  " because it has no location with at least two timepoints with non-zero effort.");
		keepMethod = false;
	}
	return std::pair(keepMethod, keepLocations);
}

void TBirpCore::_readFile(const std::string &Filename, bool hasCounts) {
	logfile().listFlushDots("Parsing ", Filename);
	// check which locations should be kept
	auto [keepMethod, keepLocations] = _keepLocations(Filename, hasCounts);
	if (!keepMethod) { return; } // chasch gd hei gah

	// open file and get relevant column indices
	coretools::TInputMaybeRcppFile file(Filename, coretools::FileType::Header, getDelim(Filename));
	auto [countIndices, covEffortIndices, covDetectionIndices, methodName] = getIndicesMethodFile(file, hasCounts);

	// assemble species names and id's
	std::vector<size_t> species_ids = {};
	if (hasCounts) {
		std::vector<std::string> speciesNames = getSpeciesOrCovariateNames(file, countIndices);
		for (const auto &name : speciesNames) { species_ids.push_back(_speciesNames.add(name)); }
	} else {
		species_ids = _simulateSpeciesNames();
	}

	// get covariate effort names and store in unique container
	std::vector<size_t> covariateEffort_ids    = fillCovNames(file, covEffortIndices, _covariateEffortNames);
	std::vector<size_t> covariateDetection_ids = fillCovNames(file, covDetectionIndices, _covariateDetectionNames);

	// add new method: name is derived from filename
	TMethods method(methodName, _data.size(), species_ids, covariateEffort_ids, covariateDetection_ids);
	_data.addMethod(method); // the index of the method is the current size of the _methods vector

	for (; !file.empty(); file.popFront()) {
		std::string location = (std::string)file.get("location");
		if (!keepLocations[location]) { continue; }                  // skip this location
		if (allEffortsAreZero(file, covEffortIndices)) { continue; } // skip timepoint: effort = 0

		size_t location_id  = _locations.add((std::string)location);
		size_t timepoint_id = _timepoints.add(file.get<TypeTime>("timepoint"));
		const auto &row     = file.front();

		// boot one vector together (stiflä) of all species counts, and one of all covariates
		std::vector<TypeCounts> countVec(species_ids.size(), 0);
		std::vector<TypeCovariateEffort> covEffortVec(covEffortIndices.size());
		std::vector<TypeCovariateDetection> covDetectionVec(covDetectionIndices.size());
		if (hasCounts) {
			for (size_t c = 0; c < countVec.size(); c++) {
				countVec[c] = coretools::str::fromStringCheck<TypeCounts>(row[countIndices[c]]);
			}
		}
		for (size_t c = 0; c < covEffortVec.size(); c++) {
			covEffortVec[c] = coretools::str::fromStringCheck<TypeCovariateEffort>(row[covEffortIndices[c]]);
		}
		for (size_t c = 0; c < covDetectionVec.size(); c++) {
			covDetectionVec[c] = coretools::str::fromStringCheck<TypeCovariateDetection>(row[covDetectionIndices[c]]);
		}
		// create timepoint
		TTimepoints pt(countVec, covEffortVec, covDetectionVec, timepoint_id);
		// add timepoint to data, either in existing object or create new one
		_data.back().add(pt, location_id); // back: goes to last element
	}
	// standardize covariates
	_data.back().standardizeCovariates();

	_data.back().fillIDVectors(_speciesNames.size());

	logfile().done();
}

void TBirpCore::_sortTimepoints() {
	// sort timepoints and corresponding locations(?) for every method/species
	// check if data (location) has more than one entry with effort > 0

	std::vector<size_t> sortingIndex = coretools::rankSort(_timepoints);
	_timepoints                      = coretools::sortContainerByRank(_timepoints, sortingIndex);
	for (auto &method : _data) { method.sorttimes(sortingIndex); }
}

void TBirpCore::_runMCMC(bool Stochastic) {
	// create model
	stattools::instances::dagBuilder().clear();
	TBirpModel model(_data, _locations, _timepoints, _speciesNames, _covariateEffortNames, _covariateDetectionNames,
	                 _outname, Stochastic, false);

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(_outname);

	// tell model to write gamma summaries to file
	model.writeSummaryGammaPosterior(_outname);
}

void TBirpCore::infer() {
	// read data
	_readForInference(true);

	if (parameters().exists("stochastic")) {
		logfile().list("Will model stochastic trends (argument 'stochastic').");
		// run MCMC (deterministic)
		logfile().startIndent("Running MCMC with deterministic model for initialization:");
		_runMCMC(false);
		logfile().endIndent();
		// infer stochastic and use parameter estimates of deterministic to initialize
		logfile().startIndent("Running MCMC with stochastic model:");
		_runMCMC(true);
		logfile().endIndent();
	} else {
		// run MCMC (deterministic)
		logfile().list("Will model deterministic trends (use 'stochastic' to change).");
		_runMCMC(false);
	}
}

//-----------------------------------
// Simulate

std::vector<std::string> TBirpCore::_getMethodspecificHeader(const TMethods &method) {
	std::vector<std::string> header = getHeader();
	// add species name to base header vector
	for (const auto &speciesID : method.speciesIDsinUniqueContainer()) {
		header.push_back("counts_" + _speciesNames[speciesID]);
	}
	// add covariate effort name to base header vector
	for (const auto &covID : method.covariateEffortIDsinUniqueContainer()) {
		header.push_back("covEffort_" + _covariateEffortNames[covID]);
	}
	// add covariate detection name to base header vector
	for (const auto &covID : method.covariateDetectionIDsinUniqueContainer()) {
		header.push_back("covDetection_" + _covariateDetectionNames[covID]);
	}
	return header;
}

void TBirpCore::_writeFile(const std::string &Name) {
	// make second file that stores names of all methodFiles generated during the simulation
	coretools::TOutputMaybeRcppFile fileList(_outname + "_" + Name + "_allFilesGenerated");
	// write file line by line
	for (auto &method : _data) {
		// get filename here.
		std::string Filename            = _outname + "_" + method.name() + "_" + Name + "_" + "counts.txt";
		// get header here
		std::vector<std::string> header = _getMethodspecificHeader(method);
		// append to fileList
		fileList << Filename << coretools::endl;
		// now create actual file with data and write
		coretools::TOutputMaybeRcppFile outfile(Filename, header);
		method.write(outfile, _locations, _timepoints);
	}
}

void TBirpCore::_writeTimePoints() {
	coretools::TOutputFile file(_outname + "_timepoints.txt", {"timepoints"});
	for (auto it : _timepoints) { file << it << coretools::endl; }
}

// todo Liam mach no ds me der abstand cha aagäh
void TBirpCore::_fillTimepointsFromCommandLine() {
	std::vector<TypeTime> timepoints;
	std::vector<TypeTime> timepointsDefault = {1.0, 2.0, 3.0};
	if (parameters().exists("timepoints")) {
		auto timepointName = parameters().get<std::string>("timepoints");
		if (coretools::str::stringIsProbablyANumber(timepointName)) {
			auto timepointNumber = coretools::str::fromString<size_t>(timepointName);
			for (size_t i = 1; i <= timepointNumber; i++) { timepoints.emplace_back((double)i); }
			if (timepoints.size() == 1) { UERROR("Please provide more than one timepoint"); }
		} else {
			parameters().fill("timepoints", timepoints);
			if (timepoints.size() == 1) { UERROR("Please provide more than one timepoint"); }
		}
	} else {
		timepoints = timepointsDefault;
	}
	logfile().list("Timepoints: ", timepoints, " (argument 'timepoints').");
	// fill into unique container
	for (const auto &timepoint : timepoints) { _timepoints.add(timepoint); }
}

void TBirpCore::_fillLocationsFromCommandLine() {
	size_t numLocations = parameters().get("numLocations", 2);
	logfile().list("Number of locations: ", numLocations, " (argument 'numLocations').");
	for (size_t l = 0; l < numLocations; ++l) { _locations.add("Location_" + coretools::str::toString(l)); }
}

std::vector<size_t> TBirpCore::_simulateSpeciesNames() {
	size_t numSpecies = parameters().get("numSpecies", 1);
	logfile().list("Number of species: ", numSpecies, " (argument 'numSpecies').");
	std::vector<size_t> species_ids(numSpecies);
	for (size_t n = 0; n < numSpecies; n++) {
		species_ids[n] = _speciesNames.add("species_" + coretools::str::toString(n + 1));
	}
	return species_ids;
}

std::vector<size_t> TBirpCore::_simulateCovEffortNames() {
	size_t numCovs = parameters().get("numCovariatesEffort", 1);
	if (numCovs == 0) { UERROR("Argument 'numCovariatesEffort' must be at least 1!"); }
	logfile().list("Number of effort covariates: ", numCovs, " (argument 'numCovariatesEffort').");

	// fill names and ids
	std::vector<size_t> covariates_ids(numCovs);
	for (size_t n = 0; n < numCovs; n++) {
		covariates_ids[n] = _covariateEffortNames.add(coretools::str::toString(n + 1));
	}
	return covariates_ids;
}

std::vector<size_t> TBirpCore::_simulateCovDetectionNames() {
	size_t numCovs = parameters().get("numCovariatesDetection", 0);
	logfile().list("Number of detection probability covariates: ", numCovs, " (argument 'numCovariatesDetection').");
	// fill names and ids
	std::vector<size_t> covariates_ids(numCovs);
	for (size_t n = 0; n < numCovs; n++) {
		covariates_ids[n] = _covariateDetectionNames.add(coretools::str::toString(n + 1));
	}
	return covariates_ids;
}

void TBirpCore::simulate() {
	// read incomplete counts file (if provided)
	if (parameters().exists("simulationFile")) {
		std::string filename = parameters().get("simulationFile");
		logfile().list("Will use simulation file '", filename, "' to simulate counts (argument 'simulationFile').");
		_readFile(filename, false);
	} else {
		// else read main dimensions from command line
		_fillTimepointsFromCommandLine();
		_fillLocationsFromCommandLine();
		size_t numMethods = parameters().get("numMethods", 1);
		logfile().list("Number of methods: ", numMethods, " (argument 'numMethods').");
		// create vector of speciesnames species_1, species_2 ... as many as there are. currently only supports one
		// species
		auto species_ids            = _simulateSpeciesNames();
		auto covariateEffort_ids    = _simulateCovEffortNames();
		auto covariateDetection_ids = _simulateCovDetectionNames();

		for (size_t m = 0; m < numMethods; ++m) {
			TMethods method("Method_" + coretools::str::toString(m + 1), _locations.size(), _timepoints.size(), m,
			                species_ids, covariateEffort_ids, covariateDetection_ids);
			_data.addMethod(method);
		}
		for (auto &m : _data) { m.fillIDVectors(_speciesNames.size()); }
	}
	_data.fillMethLocIndices(_speciesNames.size(), _locations);

	// create model
	const bool stochastic = parameters().exists("stochastic");
	if (stochastic) {
		logfile().list("Will simulate stochastic trends (argument 'stochastic').");
	} else {
		logfile().list("Will simulate deterministic trends (use 'stochastic' to change).");
	}

	TBirpModel model(_data, _locations, _timepoints, _speciesNames, _covariateEffortNames, _covariateDetectionNames,
	                 _outname, stochastic, true);

	// run simulator
	stattools::TSimulator simulator;
	simulator.simulate();

	// write Birp file
	_writeFile("simulated");
}
