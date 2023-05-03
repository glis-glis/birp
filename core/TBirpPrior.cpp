//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"
#include "TBirpCore.h"
#include "TModels.h"
#include "coretools/Storage/TDataFile.h"
#include <utility>

using namespace coretools::instances;

TBirpPrior::TBirpPrior(TData &Data, const TUniqueContainer<std::string> &Locations,
                       const TUniqueContainer<std::size_t> &TimePoints,
                       const TUniqueContainer<std::string> &CovariateNames, size_t numSpecies,
                       TBirpCommonModel &CommonModel, std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn,
                       std::unique_ptr<TBirpStochasticAddOn> &StochasticAddOn, std::string Prefix)
    : _data(Data), _locations(Locations), _timepoints(TimePoints), _covariateNames(CovariateNames),
      _numSpecies(numSpecies), _out(std::move(Prefix)) {

	// stattools: tell which ones are prior parameters
	addPriorParameter({&CommonModel.gamma, &CommonModel.timesOfChange});
	for (auto &b : CommonModel.beta) { addPriorParameter(b.get()); }
	if (_numSpecies > 1) { UERROR("Birp only supports single species at the moment.."); }
	if (NegBinomAddOn) { addPriorParameter({&NegBinomAddOn->mu, &NegBinomAddOn->logA}); }
	if (StochasticAddOn) { addPriorParameter({&StochasticAddOn->logPhi, &StochasticAddOn->sigma}); }

	// create models
	size_t speciesID = 0;
	if (StochasticAddOn) {
		_models.emplace_back(std::make_unique<TStochastic>(speciesID, Data, NegBinomAddOn, &CommonModel.gamma,
		                                                   &CommonModel.timesOfChange, &StochasticAddOn->logPhi,
		                                                   &StochasticAddOn->sigma));
	} else {
		_models.emplace_back(std::make_unique<TDeterministic>(speciesID, Data, NegBinomAddOn, &CommonModel.gamma,
		                                                      &CommonModel.timesOfChange));
	}

	for (size_t i = 0; i < Data.size(); ++i) { _data[i].addBeta(CommonModel.beta[i]); }
}

std::string TBirpPrior::name() const { return "birp"; }

//------------------------------------
// Initialization

void TBirpPrior::initializeInferred() {
	// set size of parameters etc.
	for (auto &_method : _data) { _method.initialize(_covariateNames); }
	for (size_t s = 0; s < _numSpecies; s++) { _models[s]->initialize(_data, _locations, _timepoints); }
}

void TBirpPrior::estimateInitialPriorParameters() {
	// set initial values of parameters
	for (auto &_method : _data) { _method.estimateInitialBeta(); }
	for (size_t s = 0; s < _numSpecies; ++s) { _models[s]->estimateInitial(_data, _timepoints); }

	// now initialize temporary variables (LL, JP etc)
	for (size_t s = 0; s < _numSpecies; ++s) { _models[s]->initializeTempVariables(_data); }
}

//------------------------------------
// MCMC

void TBirpPrior::updateParams() {
	for (size_t s = 0; s < _numSpecies; ++s) { _models[s]->update(_data); }
	for (auto &_method : _data) { _method.updateBeta(_data, _models); }
}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {
	// calculates likelihood
	double sum = 0.0;
	for (size_t s = 0; s < _numSpecies; ++s) { sum += _models[s]->sumCurLL(); }
	return sum;
}

//------------------------------------
// Simulations

std::vector<size_t> TBirpPrior::_readN0FromCommandLine() const {
	std::vector<size_t> N0;
	std::vector<size_t> N0default = {10};
	parameters().fillParameterIntoContainerWithDefault("N0", N0, ',', N0default);
	if (N0.size() != _numSpecies) {
		UERROR("Size of N0 (", N0.size(), ") does not match the number of species (", _numSpecies, ")!");
	}
	return N0;
}

auto TBirpPrior::_readDetectionProbFile() const {
	std::vector<std::string> methodNamesFromDetections;
	std::vector<std::string> locationNamesFromDetections;
	std::vector<coretools::Probability> detectionsProbs;
	if (parameters().parameterExists("detectionProbFile")) {
		// TODO are detection probabilities the same across all species or one file per species?
		std::string filename = parameters().getParameterFilename("detectionProbFile");

		coretools::TInputFile file(filename, coretools::TFile_Filetype::header);
		if (file.numCols() != 3 || file.header()[0] != "method" || file.header()[1] != "location" ||
		    file.header()[3] != "detectionProbs") {
			UERROR("Invalid header in file '", filename,
			       "': Expected 3 columns with names 'method', 'location' and 'detectionsProbs'!");
		}

		std::vector<std::string> tmp;
		while (file.read(tmp)) {
			methodNamesFromDetections.push_back(tmp[0]);
			locationNamesFromDetections.push_back(tmp[1]);
			detectionsProbs.push_back(coretools::str::fromString<coretools::Probability>(tmp[2]));
		}
	} else {
		// detectionProbFile is not given -> just create as many method/locations as there are and use default = 1.0
		for (auto &method : _data) {
			for (auto &loc : method) {
				methodNamesFromDetections.push_back(method.name());
				locationNamesFromDetections.push_back(_locations[loc.location_id()]);
				detectionsProbs.emplace_back(1.0);
			}
		}
	}
	return std::make_tuple(methodNamesFromDetections, locationNamesFromDetections, detectionsProbs);
}

void TBirpPrior::_simulateUnderPrior(Storage *) { // this function gets called when simulated
	// first calculate effort
	for (auto &method : _data) { method.initializeEffort(); }

	auto detectionProbs = _readDetectionProbFile();
	auto N0             = _readN0FromCommandLine();

	// do we need to simulate covariates?
	bool simulateCovariates = !parameters().parameterExists("simulationFile");

	for (size_t s = 0; s < _models.size(); s++) {
		_models[s]->simulate(_data, N0[s], detectionProbs, _locations, simulateCovariates);
	}
}

//------------------------------------
// Getters

const TUniqueContainer<std::string> &TBirpPrior::locations() const { return _locations; }

const TData &TBirpPrior::data() const { return _data; }

const TUniqueContainer<size_t> &TBirpPrior::timepoints() const { return _timepoints; }
