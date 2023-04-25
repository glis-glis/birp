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
                             const std::vector<size_t> &countIndices, const std::vector<size_t> &effortIndices) {
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
		if (2 + countIndices.size() + effortIndices.size() < 4) {
			UERROR("Header of file", Infile.name(),
			       " , should contain at least 4 columns: 'location', 'timepoint' and at least two more columns "
			       "starting with 'counts_' and 'cov_', respectively.");
		}
	} else {
		if (2 + countIndices.size() + effortIndices.size() < 3) {
			UERROR("Header of file", Infile.name(),
			       " , should contain at least 3 columns: 'method', 'location' and at least one column starting with "
			       "'cov'.");
		}
	}
}

void TBirpCore::_readForInference(bool hasCounts) {
	std::string input = parameters().getParameterFilename("data");
	bool dataIsFile =
	    std::filesystem::exists(input); // if multiple files are supplied separated by commas, then this will be false
	std::vector<std::string> filenames;
	if (dataIsFile) {
		coretools::TInputFile tryFile(input, coretools::TFile_Filetype::header);
		// first check if maybe the file is the actual counts file and not a list of files
		if (tryFile.numCols() != 1) {
			filenames.push_back(input);
		} else { // if not, then its a list of files (methods), so we fill them in
			coretools::TInputFile infile(input, coretools::TFile_Filetype::fixed); // fixed: no header
			std::vector<std::string> row;
			while (infile.read(row)) { filenames.push_back(row[0]); }
		}
	} else {
		coretools::str::fillContainerFromString(input, filenames, ',');
	}
	for (const auto &name : filenames) { _readFile(name, hasCounts); }
	_sortTimepoints();
	_data.fillMethLocIndices(_speciesNames.size());
}

std::vector<size_t> TBirpCore::getColumnIndex(const coretools::TInputFile &infile, const std::string &string) {
	std::vector<size_t> countIndices;
	for (size_t i = 0; i < infile.numCols(); i++) {
		if (coretools::str::stringStartsWith(infile.headerAt(i), string)) { countIndices.push_back(i); }
	}
	return countIndices;
}

std::vector<std::string> TBirpCore::_getSpeciesOrCovariateNames(const coretools::TInputFile &infile,
                                                                const std::vector<size_t> &Indices) {
	std::vector<std::string> names(Indices.size());
	std::string test = infile.headerAt(2);
	for (size_t s = 0; s < Indices.size(); s++) {
		names[s] = coretools::str::readAfter(infile.headerAt(Indices[s]), '_');
	}
	return names;
}

void TBirpCore::_readFile(const std::string &Filename, bool hasCounts) {
	std::string delim       = "\t";
	std::string_view ending = coretools::str::readAfterLast(Filename, '.');
	if (ending == "csv") { delim = ","; }
	coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header, delim);

	// get position of columns that contain counts: might be an empty vector in case the file has no counts
	std::vector<size_t> countIndices = getColumnIndex(infile, "counts");

	// get position of columns that contain covariates (efforts)
	std::vector<size_t> effortIndices = getColumnIndex(infile, "cov");

	_checkHeader(infile, hasCounts, countIndices, effortIndices);

	std::vector<size_t> species_ids = {};
	if (hasCounts) {
		std::vector<std::string> speciesNames = _getSpeciesOrCovariateNames(infile, countIndices);
		for (const auto &name : speciesNames) { species_ids.push_back(_speciesNames.add(name)); }
	} else {
		species_ids = _simulateSpeciesNames();
	}

	std::vector<std::string> covariateNames = _getSpeciesOrCovariateNames(infile, effortIndices);
	covariateNames.insert(covariateNames.begin(), "intercept"); // add extra intercept covariate
	std::vector<size_t> covariate_ids(covariateNames.size());
	for (size_t c = 0; c < covariateNames.size(); ++c) { covariate_ids[c] = _covariateNames.add(covariateNames[c]); }

	// add new method: name is derived from filename
	std::string_view methodNameFull = Filename;
	if (coretools::str::stringContains(Filename, '/')) {
		methodNameFull = coretools::str::readAfterLast(Filename, '/');
	}
	std::string_view methodNameCropped = coretools::str::readAfter(methodNameFull, '_');
	std::string methodName             = (std::string)coretools::str::readBeforeLast(methodNameCropped, '.');
	coretools::str::eraseAllOccurences(methodName, "_simulatedCounts");
	TMethods method(methodName, _data.size(), species_ids, covariate_ids);
	_data.addMethod(method); // the index of the method is the current size of the _methods vector

	std::vector<std::string> row;
	while (infile.read(row)) {
		size_t location_id  = _locations.add(row[0]);
		size_t timepoint_id = _timepoints.add(coretools::str::fromStringCheck<size_t>(row[1]));

		// boot one vector together (stiflä) of all species counts, and one of all covariates
		std::vector<TypeCounts> countVec(species_ids.size(), 0);
		std::vector<TypeCovariate> covVec(effortIndices.size());
		if (hasCounts) {
			for (size_t c = 0; c < countVec.size(); c++) { countVec[c] = (TypeCounts)row[countIndices[c]]; }
		}
		for (size_t c = 0; c < covVec.size(); c++) { covVec[c] = (TypeCovariate)row[effortIndices[c]]; }
		// create timepoint
		TTimepoints pt(countVec, covVec, timepoint_id);
		// add timepoint to data, either in existing object or create new one
		_data.back().add(pt, location_id); // back: goes to last element
	}

	// standardize covariates
	_data.back().standardizeCovariatesAndAddIntercept();

	_data.back().fillIDVectors(_speciesNames.size());
}

void TBirpCore::_sortTimepoints() {
	// sort timepoints and corresponding locations(?) for every method/species
	// check if data (location) has more than one entry with effort > 0

	std::vector<size_t> sortingIndex = coretools::rankSort(_timepoints);
	_timepoints                      = coretools::sortContainerByRank(_timepoints, sortingIndex);
	for (auto &method : _data) {
		method.sorttimes(sortingIndex);
	}
}

void TBirpCore::infer() {
	// read data
	_readForInference(true);

	// create model
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_data, _locations, _timepoints, _speciesNames, _covariateNames, _outname, dagBuilder);

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(_outname, &dagBuilder);
}

//-----------------------------------
// Simulate

std::vector<std::string> TBirpCore::_getMethodspecificHeader(const TMethods &method) {
	std::vector<std::string> header = _getHeader();
	// add species name to base header vector
	for (const auto &speciesID : method.speciesIDsinUniqueContainer()) {
		header.push_back("counts_" + _speciesNames[speciesID]);
	}
	// add covariate name to base header vector
	for (const auto &covID : method.covariateIDsinUniqueContainer()) {
		if (_covariateNames[covID] == "intercept") { continue; } // do not write intercept to file
		header.push_back("cov_" + _covariateNames[covID]);
	}
	return header;
}

void TBirpCore::_writeFile() {
	// make second file that stores names of all methodFiles generated during the simulation
	coretools::TOutputFile fileList(_outname + "_allFilesGenerated");
	// write file line by line
	for (auto &method : _data) {
		// get filename here.
		std::string Filename            = _outname + "_" + method.name() + "_simulatedCounts.txt";
		// get header here
		std::vector<std::string> header = _getMethodspecificHeader(method);
		// append to fileList
		fileList << _outname + "_" + method.name() + "_simulatedCounts.txt"  << coretools::endl;
		// now create actual file with data and write
		coretools::TOutputFile outfile(Filename, header);
		method.write(outfile, _locations, _timepoints);
	}
}

// todo lol Madi lueg ds ah
void TBirpCore::_fillTimepointsFromCommandLine() {
	std::vector<size_t> timepoints;
    std::vector<size_t> timepointsDefault = {1, 2, 3};
    if (parameters().parameterExists("timepoints")){
        auto timepointName = parameters().getParameter<std::string>("timepoints");
        if (coretools::str::stringIsProbablyANumber(timepointName)) {
            size_t timepointNumber = coretools::str::fromString<size_t>(timepointName);
            for(size_t i = 1; i <= timepointNumber; i++){
                timepoints.push_back(i);
            }
        } else {
            parameters().fillParameterIntoContainer("timepoints", timepoints, ',');
        }
    } else {
        parameters().fillParameterIntoContainerWithDefault("timepoints", timepoints, ',', timepointsDefault);
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

std::vector<size_t> TBirpCore::_simulateCovNames() {
	size_t numCovs = parameters().getParameterWithDefault("numCovariates", 1);
	if (numCovs == 0) { UERROR("Argument 'numCovariates' must be at least 1!"); }
	// fill names and ids and add intercept
	std::vector<size_t> covariates_ids(numCovs + 1);
	covariates_ids[0] = _covariateNames.add("intercept");
	for (size_t n = 1; n < numCovs + 1; n++) { covariates_ids[n] = _covariateNames.add(coretools::str::toString(n)); }
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
		size_t numMethods  = parameters().getParameterWithDefault("numMethods", 1);
		// create vector of speciesnames species_1, species_2 ... as many as there are. currently only supports one
		// species
		auto species_ids   = _simulateSpeciesNames();
		auto covariate_ids = _simulateCovNames();

		for (size_t m = 0; m < numMethods; ++m) {
			TMethods method("Method_" + coretools::str::toString(m + 1), _locations.size(), _timepoints.size(), m,
			                species_ids, covariate_ids);
			_data.addMethod(method);
		}
		for (auto &m : _data) { m.fillIDVectors(_speciesNames.size()); }
	}
	_data.fillMethLocIndices(_speciesNames.size());

	// create model
	stattools::TDAGBuilder dagBuilder;
	TBirpModel model(_data, _locations, _timepoints, _speciesNames, _covariateNames, _outname, dagBuilder);

	// run simulator
	stattools::TSimulator simulator;
	simulator.simulate(_outname, &dagBuilder); // goes to ... DAG ... TBirpPrior

	// write Birp file
	_writeFile();
}
