/*
 * TBirpCore.cpp
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#include "TBirpCore.h"
#include "TData.h"
#include "coretools/Types/strongTypes.h"
#include "stattools/MCMC/TMCMC.h"
#include "stattools/Priors/TPriorDirichlet.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorNormal.h"

using namespace coretools::instances;
using namespace stattools;
using namespace stattools::prior;

//--------------------------------------
// TBirpCore
//--------------------------------------

TBirpCore::TBirpCore() {
	using namespace coretools::instances;
	if (parameters().parameterExists("out")) {
		_outname = parameters().getParameter<std::string>("out");
		logfile().list("Will write files with prefix '", _outname, "'. (parameter 'out')");
	} else {
		_outname = "birp";
		logfile().list("Will write files with default prefix '", _outname,
		               "'. (You can change this using the parameter 'out')");
	}
}

//-----------------------------------
// Infer

bool TBirpCore::_checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index,
                                      const std::string &Title) const {
	return Infile.headerAt(index) == Title;
}

std::vector<std::string> TBirpCore::_getHeader() const { return {"location", "timepoint"}; }

void TBirpCore::_checkHeader(const coretools::TInputFile &Infile, bool hasCounts,
                             const std::vector<size_t> &countIndices,
                             const std::vector<size_t> &CovEffortIndices) const {
	std::vector<std::string> headerVec = _getHeader();
	// check whether first two columns are correct
	for (size_t i = 0; i < headerVec.size(); i++) {
		if (!_checkHeaderPlacement(Infile, i, headerVec[i])) {
			UERROR("Error in header of file", Infile.name(), " , please make sure header at column ", i + 1, " is ",
			       headerVec[i], ". Found: '", Infile.headerAt(i), "'");
		}
	}
	// check if number of columns is at least 4 (with counts) or 3 (without counts)
	if (hasCounts) {
		if (headerVec.size() + countIndices.size() + CovEffortIndices.size() < 4) {
			UERROR("Header of file", Infile.name(),
			       " , should contain at least 4 columns: 'location', 'timepoint' and at least two more columns "
			       "starting with 'counts_' and 'covEffort_', respectively.");
		}
	} else {
		if (headerVec.size() + countIndices.size() + CovEffortIndices.size() < 3) {
			UERROR(
			    "Header of file", Infile.name(),
			    " , should contain at least 3 columns: 'location', 'timepoint' and at least one column starting with "
			    "'covEffort_'.");
		}
	}
}

std::string TBirpCore::_getDelim(const std::string &Filename) const {
	std::string delim       = "\t";
	std::string_view ending = coretools::str::readAfterLast(Filename, '.');
	if (ending == "csv") { delim = ","; }
	return delim;
}

void TBirpCore::_readForInference(bool hasCounts) {
	std::string input = parameters().getParameterFilename("data");
	bool dataIsFile =
	    std::filesystem::exists(input); // if multiple files are supplied separated by commas, then this will be false
	std::vector<std::string> filenames;
	if (dataIsFile) {
		coretools::TInputFile tryFile(input, coretools::TFile_Filetype::header, _getDelim(input));
		// first check if maybe the file is the actual counts file and not a list of files
		if (tryFile.numCols() != 1) {
			filenames.push_back(input);
		} else { // if not, then its a list of files (methods), so we fill them in
			coretools::TInputFile infile(input, coretools::TFile_Filetype::fixed, _getDelim(input)); // fixed: no header
			std::vector<std::string> row;
			while (infile.read(row)) { filenames.push_back(row[0]); }
		}
	} else {
		coretools::str::fillContainerFromString(input, filenames, ',');
	}
	for (const auto &name : filenames) { _readFile(name, hasCounts); }
	_sortTimepoints();
	_data.fillMethLocIndices(_speciesNames.size());

	// write filtered counts
	_writeFile("filtered");
}

std::vector<size_t> TBirpCore::getColumnIndex(const coretools::TInputFile &infile, const std::string &string) const {
	std::vector<size_t> countIndices;
	for (size_t i = 0; i < infile.numCols(); i++) {
		if (coretools::str::stringStartsWith(infile.headerAt(i), string)) { countIndices.push_back(i); }
	}
	return countIndices;
}

std::vector<std::string> getSpeciesOrCovariateNames(const coretools::TInputFile &infile,
                                                    const std::vector<size_t> &Indices) {
	std::vector<std::string> names(Indices.size());
	std::string test = infile.headerAt(2);
	for (size_t s = 0; s < Indices.size(); s++) {
		names[s] = coretools::str::readAfter(infile.headerAt(Indices[s]), '_');
	}
	return names;
}

std::vector<size_t> fillCovNames(coretools::TInputFile &infile, const std::vector<size_t> &Indices,
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
	std::string_view methodNameCropped = coretools::str::readAfter(methodNameFull, '_');
	std::string methodName             = (std::string)coretools::str::readBeforeLast(methodNameCropped, '.');
	coretools::str::eraseAllOccurences(methodName, "_simulatedCounts");
	return methodName;
}

auto TBirpCore::_getIndicesMethodFile(const coretools::TInputFile &Infile, bool hasCounts) const {
	// get position of columns that contain counts: might be an empty vector in case the file has no counts
	std::vector<size_t> countIndices = getColumnIndex(Infile, "counts");

	// get position of columns that contain covariates (efforts)
	std::vector<size_t> covEffortIndices    = getColumnIndex(Infile, "covEffort");
	std::vector<size_t> covDetectionIndices = getColumnIndex(Infile, "covDetection");

	_checkHeader(Infile, hasCounts, countIndices, covEffortIndices);

	auto methodName = getMethodName(Infile.name());
	return std::make_tuple(countIndices, covEffortIndices, covDetectionIndices, methodName);
}

bool TBirpCore::_allEffortsAreZero(const std::vector<std::string> &Row, const std::vector<size_t> &CovEffortIndices,
                                   const std::string &Filename) const {
	std::set<TypeCovariateEffort> covEffort;
	for (auto c : CovEffortIndices) { covEffort.insert((TypeCovariateEffort)Row[c]); }

	if (covEffort.empty()) {
		UERROR("Invalid file ", Filename, ": Need at least one column starting with 'covEffort_'.");
	}

	return covEffort.size() == 1 && *covEffort.begin() == 0.0; // returns true if all covariates of effort are zero
}

auto TBirpCore::_keepLocations(const std::string &Filename, bool hasCounts) const {
	// open file and get relevant column indices
	coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header, _getDelim(Filename));
	auto [countIndices, covEffortIndices, covDetectionIndices, methodName] = _getIndicesMethodFile(infile, hasCounts);

	// map: per location, how many timepoints with non-zero effort are there?
	std::map<std::string, size_t> numTimePointsPerLocation;
	std::vector<std::string> row;
	while (infile.read(row)) {
		// read location and timepoint
		const std::string location  = row[0];
		const std::string timepoint = row[1];

		// read effort(s)
		if (_allEffortsAreZero(row, covEffortIndices, Filename)) {
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
	// check which locations should be kept
	auto [keepMethod, keepLocations] = _keepLocations(Filename, hasCounts);
	if (!keepMethod) { return; } // chasch gd hei gah

	// open file and get relevant column indices
	coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header, _getDelim(Filename));
	auto [countIndices, covEffortIndices, covDetectionIndices, methodName] = _getIndicesMethodFile(infile, hasCounts);

	// assemble species names and id's
	std::vector<size_t> species_ids = {};
	if (hasCounts) {
		std::vector<std::string> speciesNames = getSpeciesOrCovariateNames(infile, countIndices);
		for (const auto &name : speciesNames) { species_ids.push_back(_speciesNames.add(name)); }
	} else {
		species_ids = _simulateSpeciesNames();
	}

	// get covariate effort names and store in unique container
	std::vector<size_t> covariateEffort_ids    = fillCovNames(infile, covEffortIndices, _covariateEffortNames);
	std::vector<size_t> covariateDetection_ids = fillCovNames(infile, covDetectionIndices, _covariateDetectionNames);

	// add new method: name is derived from filename
	TMethods method(methodName, _data.size(), species_ids, covariateEffort_ids, covariateDetection_ids);
	_data.addMethod(method); // the index of the method is the current size of the _methods vector

	std::vector<std::string> row;
	while (infile.read(row)) {
		const std::string location = row[0];
		if (!keepLocations[location]) { continue; }                            // skip this location
		if (_allEffortsAreZero(row, covEffortIndices, Filename)) { continue; } // skip this timepoint: effort = 0

		size_t location_id  = _locations.add(location);
		size_t timepoint_id = _timepoints.add(coretools::str::fromStringCheck<size_t>(row[1]));

		// boot one vector together (stiflä) of all species counts, and one of all covariates
		std::vector<TypeCounts> countVec(species_ids.size(), 0);
		std::vector<TypeCovariateEffort> covEffortVec(covEffortIndices.size());
		std::vector<TypeCovariateDetection> covDetectionVec(covDetectionIndices.size());
		if (hasCounts) {
			for (size_t c = 0; c < countVec.size(); c++) { countVec[c] = (TypeCounts)row[countIndices[c]]; }
		}
		for (size_t c = 0; c < covEffortVec.size(); c++) {
			covEffortVec[c] = (TypeCovariateEffort)row[covEffortIndices[c]];
		}
		for (size_t c = 0; c < covDetectionVec.size(); c++) {
			covDetectionVec[c] = (TypeCovariateDetection)row[covDetectionIndices[c]];
		}
		// create timepoint
		TTimepoints pt(countVec, covEffortVec, covDetectionVec, timepoint_id);
		// add timepoint to data, either in existing object or create new one
		_data.back().add(pt, location_id); // back: goes to last element
	}

	// standardize covariates
	_data.back().standardizeCovariates();

	_data.back().fillIDVectors(_speciesNames.size());
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
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_data, _locations, _timepoints, _speciesNames, _covariateEffortNames, _covariateDetectionNames,
	                 _outname, dagBuilder, Stochastic);

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(_outname, &dagBuilder);

	// tell model to write gamma summaries to file
	model.writeSummaryGammaPosterior(_outname);
}

void TBirpCore::infer() {
	// read data
	_readForInference(true);

	if (parameters().parameterExists("stochastic")) {
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
		_runMCMC(false);
	}
}

//-----------------------------------
// Simulate

std::vector<std::string> TBirpCore::_getMethodspecificHeader(const TMethods &method) {
	std::vector<std::string> header = _getHeader();
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
	coretools::TOutputFile fileList(_outname + "_" + Name + "_allFilesGenerated");
	// write file line by line
	for (auto &method : _data) {
		// get filename here.
		std::string Filename            = _outname + "_" + method.name() + "_" + Name + "_" + "counts.txt";
		// get header here
		std::vector<std::string> header = _getMethodspecificHeader(method);
		// append to fileList
		fileList << Filename << coretools::endl;
		// now create actual file with data and write
		coretools::TOutputFile outfile(Filename, header);
		method.write(outfile, _locations, _timepoints);
	}
}

// todo Liam mach no ds me der abstand cha aagäh
void TBirpCore::_fillTimepointsFromCommandLine() {
	std::vector<size_t> timepoints;
	std::vector<size_t> timepointsDefault = {1, 2, 3};
	if (parameters().parameterExists("timepoints")) {
		auto timepointName = parameters().getParameter<std::string>("timepoints");
		if (coretools::str::stringIsProbablyANumber(timepointName)) {
			size_t timepointNumber = coretools::str::fromString<size_t>(timepointName);
			for (size_t i = 1; i <= timepointNumber; i++) { timepoints.push_back(i); }
			if (timepoints.size() == 1) { UERROR("Please provide more than one timepoint"); }
		} else {
			parameters().fillParameterIntoContainer("timepoints", timepoints, ',');
			if (timepoints.size() == 1) { UERROR("Please provide more than one timepoint"); }
		}
	} else {
		timepoints = timepointsDefault;
	}
	// fill into unique container
	for (const auto &timepoint : timepoints) { _timepoints.add(timepoint); }
}

void TBirpCore::_fillLocationsFromCommandLine() {
	size_t numLocations = parameters().getParameterWithDefault("numLocations", 2);
	for (size_t l = 0; l < numLocations; ++l) { _locations.add("Location_" + coretools::str::toString(l)); }
}

std::vector<size_t> TBirpCore::_simulateSpeciesNames() {
	size_t numSpecies = parameters().getParameterWithDefault("numSpecies", 1);
	std::vector<size_t> species_ids(numSpecies);
	for (size_t n = 0; n < numSpecies; n++) {
		species_ids[n] = _speciesNames.add("species_" + coretools::str::toString(n + 1));
	}
	return species_ids;
}

std::vector<size_t> TBirpCore::_simulateCovEffortNames() {
	size_t numCovs = parameters().getParameterWithDefault("numCovariatesEffort", 1);
	if (numCovs == 0) { UERROR("Argument 'numCovariatesEffort' must be at least 1!"); }
	// fill names and ids
	std::vector<size_t> covariates_ids(numCovs);
	for (size_t n = 0; n < numCovs; n++) {
		covariates_ids[n] = _covariateEffortNames.add(coretools::str::toString(n + 1));
	}
	return covariates_ids;
}

std::vector<size_t> TBirpCore::_simulateCovDetectionNames() {
	size_t numCovs = parameters().getParameterWithDefault("numCovariatesDetection", 0);
	// fill names and ids
	std::vector<size_t> covariates_ids(numCovs);
	for (size_t n = 0; n < numCovs; n++) {
		covariates_ids[n] = _covariateDetectionNames.add(coretools::str::toString(n + 1));
	}
	return covariates_ids;
}

void TBirpCore::simulate() {
	// read incomplete counts file (if provided)
	if (parameters().parameterExists("simulationFile")) {
		std::string filename = parameters().getParameterFilename("simulationFile");
		_readFile(filename, false);
	} else {
		// else read main dimensions from command line
		_fillTimepointsFromCommandLine();
		_fillLocationsFromCommandLine();
		size_t numMethods           = parameters().getParameterWithDefault("numMethods", 1);
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
	_data.fillMethLocIndices(_speciesNames.size());

	// create model
	const bool stochastic = parameters().parameterExists("stochastic");
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_data, _locations, _timepoints, _speciesNames, _covariateEffortNames, _covariateDetectionNames,
	                 _outname, dagBuilder, stochastic);

	// run simulator
	stattools::TSimulator simulator;
	simulator.simulate(_outname, &dagBuilder); // goes to ... DAG ... TBirpPrior

	// write Birp file
	_writeFile("simulated");
}
