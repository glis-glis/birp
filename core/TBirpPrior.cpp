//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"
#include "TEpochModel.h"
#include "coretools/Storage/TDataFile.h"
#include <utility>

using namespace coretools::instances;

TBirpPrior::TBirpPrior(std::vector<TMethods> *Methods, TUniqueContainer<std::string> *Locations,
                       TUniqueContainer<std::size_t> *TimePoints, TUniqueContainer<std::string> *CovariateNames,
                       size_t numSpecies, stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                       stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
                       std::vector<std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>>> &Beta, std::string Prefix)
    : _methods(Methods), _locations(Locations), _timepoints(TimePoints), _covariateNames(CovariateNames),
      _numSpecies(numSpecies), _out(std::move(Prefix)) {
	// constructor for Poisson model

	addPriorParameter({Gamma, TimesOfChange});
	for (auto &b : Beta) { addPriorParameter(b.get()); }
	if (_numSpecies > 1) { UERROR("Birp only supports single species at the moment.."); }
	_epochModel.emplace_back(std::make_unique<TPoissonEpochModel>(Gamma, TimesOfChange));
	for (size_t i = 0; i < _methods->size(); ++i) { (*_methods)[i].addBeta(Beta[i]); }
}

TBirpPrior::TBirpPrior(std::vector<TMethods> *Methods, TUniqueContainer<std::string> *Locations,
                       TUniqueContainer<std::size_t> *TimePoints, TUniqueContainer<std::string> *CovariateNames,
                       size_t numSpecies, stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                       stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
                       stattools::TParameterTyped<TypeMu, 1> *Mu, stattools::TParameterTyped<TypeLogA, 1> *LogA,
                       std::vector<std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>>> &Beta, std::string Prefix)
    : _methods(Methods), _locations(Locations), _timepoints(TimePoints), _covariateNames(CovariateNames),
      _numSpecies(numSpecies), _out(std::move(Prefix)) {
	// constructor for Negative Binomial model

	addPriorParameter({Gamma, TimesOfChange, Mu, LogA});
	for (auto &b : Beta) { addPriorParameter(b.get()); }
	if (_numSpecies > 1) { UERROR("Birp only supports single species at the moment.."); }
	_epochModel.emplace_back(std::make_unique<TNegBinEpochModel>(Gamma, TimesOfChange, Mu, LogA));
	for (size_t i = 0; i < _methods->size(); ++i) { (*_methods)[i].addBeta(Beta[i]); }
}

std::string TBirpPrior::name() const { return "birp"; }

//------------------------------------
// Initialization

void TBirpPrior::initializeInferred() {
	for (size_t s = 0; s < _numSpecies; s++) {
		_epochModel[s]->initializeInferred(*_timepoints, *_locations, *_methods);
	}
	for (auto &_method : *_methods) { _method.initializeInferred(*_covariateNames); }
}

void TBirpPrior::estimateInitialPriorParameters() {
	for (auto &_method : *_methods) { _method.estimateInitialBeta(); }
	for (size_t s = 0; s < _numSpecies; ++s) { _epochModel[s]->estimateInitialPriorParameters(*_methods); }
}

//------------------------------------
// MCMC

void TBirpPrior::updateParams() {
	for (size_t s = 0; s < _numSpecies; ++s) { _epochModel[s]->update(*_methods); }
	for (auto &_method : *_methods) { _method.updateBeta(_epochModel); }
}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {
	double sum = 0.0;
	for (size_t s = 0; s < _numSpecies; ++s) { sum += _epochModel[s]->calculateLL(*_methods); }
	return sum;
}

//------------------------------------
// Simulations

std::vector<size_t> TBirpPrior::_readN0FromCommandLine() const {
	std::vector<size_t> N0;
	std::vector<size_t> N0default = {100};
	parameters().fillParameterIntoContainerWithDefault("N0", N0, ',', N0default);
	return N0;
}

coretools::TMultiDimensionalStorage<double, 2> TBirpPrior::_readDetectionProbFile() const {
	if (parameters().parameterExists("detectionProbFile")) {
		// TODO are detection probabilities the same across all species or one file per species?
		std::string filename = parameters().getParameterFilename("detectionProbFile");

		auto locationNames = std::make_shared<coretools::TNamesStrings>();
		auto methodNames   = std::make_shared<coretools::TNamesStrings>();

		// read file
		coretools::TDataReaderBase<double, 2> reader;
		size_t guessNumRows = 50; // doesn't matter, not too many anyways
		coretools::TMultiDimensionalStorage<double, 2> detectionProbs;
		reader.read(filename, &detectionProbs, {locationNames, methodNames}, {'\n', '\t'}, {true, true}, {0, 0},
		            {false, false}, guessNumRows, coretools::colNames_multiline);
		return detectionProbs;
	} else {
		coretools::TMultiDimensionalStorage<double, 2> storage({_locations->size(), _methods->size()});
		storage.fill(1.0);
		auto locationNames = std::make_shared<coretools::TNamesStrings>();
		auto methodNames   = std::make_shared<coretools::TNamesStrings>();
		for (const auto &location : *_locations) { locationNames->addName({location}); }
		for (const auto &method : *_methods) { methodNames->addName({method.name()}); }
		storage.setDimensionNames({locationNames, methodNames});
		return storage;
	}
}

void TBirpPrior::_simulateUnderPrior(Storage *) { // this function gets called when simulated
	// first calculate effort
	for (auto &method : *_methods) { method.initializeEffort(); }

	if (parameters().parameterExists("simulationFile")) { // already read simulation parameters from file?
		_simulateFromFile();
	} else {
		_simulateFromCommandLine();
	}
}

void TBirpPrior::_simulateFromFile() {
	auto detectionProbs = _readDetectionProbFile();
	auto N0             = _readN0FromCommandLine();
	for (size_t s = 0; s < _epochModel.size(); s++) {
		_epochModel[s]->simulate(*_methods, N0[s], detectionProbs, *_locations);
	}
}

void TBirpPrior::_simulateFromCommandLine() {
	for (size_t s = 0; s < _epochModel.size(); s++) { _epochModel[s]->simulate(*_methods); }
}

//------------------------------------
// Getters

const TUniqueContainer<std::string> &TBirpPrior::locations() const { return *_locations; }

const std::vector<TMethods> &TBirpPrior::methods() const { return *_methods; }

const TUniqueContainer<size_t> &TBirpPrior::timepoints() const { return *_timepoints; }
