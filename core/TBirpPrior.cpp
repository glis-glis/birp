//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"

#include "TEpochModel.h"
#include "coretools/Storage/TDataFile.h"
#include <utility>

using namespace coretools::instances;

TBirpPrior::TBirpPrior(std::vector<TSpecies> *Species, TUniqueContainer<std::string> *Locations,
                       TUniqueContainer<std::string> *Methods, TUniqueContainer<std::size_t> *TimePoints,
                       stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                       stattools::TParameterTyped<TypeTime, 1> *TimesOfChange, std::string Prefix)
    : _species(Species), _locations(Locations), _methods(Methods), _timepoints(TimePoints), _out(std::move(Prefix)) {
	// constructor for Poisson model

	addPriorParameter({Gamma, TimesOfChange});
	_epochModel.emplace_back(Gamma, TimesOfChange);
}

std::string TBirpPrior::name() const { return "birp"; }

//------------------------------------
// Initialization

void TBirpPrior::initializeInferred() {
	for (size_t s = 0; s < _epochModel.size(); ++s) { _epochModel[s].initializeInferred(*_timepoints); }
}

void TBirpPrior::estimateInitialPriorParameters() {
	for (size_t s = 0; s < _species->size(); ++s) { _epochModel[s].estimateInitialPriorParameters((*_species)[s]); }
}

//------------------------------------
// MCMC

void TBirpPrior::updateParams() {
	for (size_t s = 0; s < _species->size(); ++s) { _epochModel[s].updateGamma((*_species)[s]); }
}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {
	// TODO: is this correct for multi-species case?
	double sum = 0.0;
	for (size_t s = 0; s < _species->size(); ++s) { sum += _epochModel[s].calculateLL((*_species)[s]); }
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
		for (const auto &method : *_methods) { methodNames->addName({method}); }
		storage.setDimensionNames({locationNames, methodNames});
		return storage;
	}
}

void TBirpPrior::_simulateUnderPrior(Storage *) {
	if (parameters().parameterExists("simulationFile")) { // already read simulation parameters from file?
		_simulateFromFile();
	} else {
		_simulateFromCommandLine();
	}
}

void TBirpPrior::_simulateFromFile() {
	auto detectionProbs = _readDetectionProbFile();
	auto N0             = _readN0FromCommandLine();
	for (size_t i = 0; i < _epochModel.size(); i++) {
		_epochModel[i].simulate((*_species)[i], N0[i], detectionProbs, *_locations, *_methods);
	}
}

void TBirpPrior::_simulateFromCommandLine() {
	for (size_t i = 0; i < _epochModel.size(); i++) { _epochModel[i].simulate((*_species)[i]); }
}

//------------------------------------
// Getters

const TUniqueContainer<std::string> &TBirpPrior::locations() const { return *_locations; }

const TUniqueContainer<std::string> &TBirpPrior::methods() const { return *_methods; }

const TUniqueContainer<size_t> &TBirpPrior::timepoints() const { return *_timepoints; }

const std::vector<TSpecies> &TBirpPrior::species() const { return *_species; }
