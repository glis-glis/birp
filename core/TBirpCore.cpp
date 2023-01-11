/*
 * TBirpCore.cpp
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#include "TBirpCore.h"
#include "coretools/Types/strongTypes.h"
#include "stattools/MCMC/TMCMC.h"

using namespace coretools::instances;
using namespace stattools;
using namespace stattools::prior;

TBirpModel::TBirpModel(std::vector<TSpecies> *Species, TUniqueContainer<std::string> *Locations,
                       TUniqueContainer<std::string> *Methods, TUniqueContainer<std::size_t> *TimePoints,
                       const std::string &Filename, stattools::TDAGBuilder &DAGBuilder, bool NegativeBinomial)
    // TODO: change to improper prior!
    : _gamma("gamma", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeGamma, 1>>(),
             {Filename}),
      _timesOfChange("timesOfChange",
                     std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeTime, 1>>(),
                     {Filename}),
      _mu("mu", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeMu, 2>>(), {Filename}) {

	// currently only allow for single species
	if (Species->size() != 1) { UERROR("Currently only supports single species!"); }

	// create birp prior (Poisson or Negative Binomial)
	std::shared_ptr<TBirpPrior> birp;
	if (NegativeBinomial) {
		birp = std::make_shared<TBirpPrior>(Species, Locations, Methods, TimePoints, &_gamma, &_timesOfChange, &_mu,
		                                    Filename);
	} else {
		birp =
		    std::make_shared<TBirpPrior>(Species, Locations, Methods, TimePoints, &_gamma, &_timesOfChange, Filename);
	}

	// create fake storage: data are stored in TSpecies, but stattools requires TMultiDimensionalStorage
	coretools::TMultiDimensionalStorage<TypeCounts, 2> fakeData;
	_data = std::make_unique<stattools::TObservationTyped<TypeCounts, 2>>("data", birp, fakeData,
	                                                                      stattools::TObservationDefinition(Filename));

	// add to DAG
	DAGBuilder.addToDAG({&_gamma, &_timesOfChange});
	if (NegativeBinomial) { DAGBuilder.addToDAG(&_mu); }
	DAGBuilder.addToDAG(_data.get());
	DAGBuilder.buildDAG();
};

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

std::vector<std::string> TBirpCore::_getHeader() const {
	return {"species", "location", "method", "time", "effort", "counts"};
}

void TBirpCore::_checkHeader(const coretools::TInputFile &Infile) {
	std::vector<std::string> headerVec = _getHeader();
	for (size_t i = 0; i < headerVec.size(); i++) {
		if (!_checkHeaderPlacement(Infile, i, headerVec[i])) {
			UERROR("Error in header of file", Infile.name(), " , please make sure header at column ", i + 1, " is ",
			       headerVec[i], ". Found: '", Infile.headerAt(i), "'");
		}
	}
}

void TBirpCore::_readForInference() {
	std::string filenameCounts = parameters().getParameterFilename("counts");
	_readFile(filenameCounts, true);
	_sortTimepoints();
}

void TBirpCore::_readFile(const std::string &Filename, bool hasCounts) {
	std::string delim       = "\t";
	std::string_view ending = coretools::str::readAfterLast(Filename, '.');
	if (ending == "csv") { delim = ","; }
	coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header, delim);
	if (hasCounts && infile.numCols() != 6) {
		UERROR("Expected file ", infile.name(), " to have 6 columns, but found ", infile.numCols(), ".");
	} else if (!hasCounts && infile.numCols() != 5) {
		UERROR("Expected file ", infile.name(), " to have 5 columns, but found ", infile.numCols(), ".");
	}
	_checkHeader(infile);

	std::vector<std::string> vec;
	while (infile.read(vec)) {
		size_t locationIndex  = _locations.add(vec[1]);
		size_t methodIndex    = _methods.add(vec[2]);
		size_t timepointIndex = _timepoints.add(coretools::str::fromStringCheck<size_t>(vec[3]));
		size_t counts         = 0;
		if (hasCounts) { counts = coretools::str::fromStringCheck<size_t>(vec[5]); }
		TDataPoint pt(timepointIndex, counts, coretools::str::fromStringCheck<double>(vec[4]));
		// loop over all species and check if it already exists
		bool speciesExists = false;
		for (auto &species : _species) {
			if (species.name() == vec[0]) { // did find it
				speciesExists = true;
				species.add(pt, locationIndex, methodIndex);
				break;
			}
		}
		if (!speciesExists) { _species.emplace_back(vec[0], pt, locationIndex, methodIndex); }
	}
}

void TBirpCore::_sortTimepoints() {
	// sort timepoints and corresponding dataseries for every species
	// check if data has more than one entry with effort > 0
	for (auto &species : _species) { species.checkData(_locations, _methods); }
	std::vector<size_t> sortingIndex = coretools::rankSort(_timepoints);
	_timepoints                      = coretools::sortContainerByRank(_timepoints, sortingIndex);
	for (auto &species : _species) { species.sorttimes(sortingIndex); }
}

void TBirpCore::infer() {
	// read data
	_readForInference();

	// create model
	stattools::TDAGBuilder dagBuilder;
	bool negBinom = parameters().parameterExists("negativeBinomial");
	TBirpModel model(&_species, &_locations, &_methods, &_timepoints, _outname, dagBuilder, negBinom);

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(_outname, &dagBuilder);
}

//-----------------------------------
// Simulate

void TBirpCore::_writeFile(const std::string &Filename) {
	auto header = _getHeader();
	coretools::TOutputFile outfile(Filename, header);
	// write file line by line
	for (auto &species : _species) { species.write(outfile, _locations, _methods, _timepoints); }
}

void TBirpCore::_fillTimepointsFromCommandLine() {
	std::vector<size_t> timepoints;
	std::vector<size_t> timepointsDefault = {1, 2, 3};
	parameters().fillParameterIntoContainerWithDefault("timepoints", timepoints, ',', timepointsDefault);

	// fill into unique container
	for (const auto &timepoint : timepoints) { _timepoints.add(timepoint); }
}

void TBirpCore::_fillLocationsFromCommandLine() {
	size_t numLocations = parameters().getParameterWithDefault("numLocations", 2);
	for (size_t l = 0; l < numLocations; ++l) { _locations.add("Location_" + coretools::str::toString(l)); }
}

void TBirpCore::_fillMethodsFromCommandLine() {
	size_t numMethods = parameters().getParameterWithDefault("numMethods", 2);
	for (size_t m = 0; m < numMethods; ++m) { _methods.add("Method_" + coretools::str::toString(m)); }
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
		_fillMethodsFromCommandLine();
		_species.emplace_back("Seal", _locations.size(), _methods.size(), _timepoints.size());
	}

	// create model
	stattools::TDAGBuilder dagBuilder;
	bool negBinom = parameters().parameterExists("negativeBinomial");
	TBirpModel model(&_species, &_locations, &_methods, &_timepoints, _outname, dagBuilder, negBinom);

	// run simulator
	stattools::TSimulator simulator;
	simulator.simulate(_outname, &dagBuilder);

	// write Birp file
	_writeFile(_outname + "_simulatedCounts.txt");
}
