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
	std::vector<size_t> N0default = {1000};
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
		for (auto &method : *_methods) {
			for (auto &loc : method) {
				methodNamesFromDetections.push_back(method.name());
				locationNamesFromDetections.push_back((*_locations)[loc.location_id()]);
				detectionsProbs.emplace_back(1.0);
			}
		}
	}
	return std::make_tuple(methodNamesFromDetections, locationNamesFromDetections, detectionsProbs);
}

void TBirpPrior::_simulateUnderPrior(Storage *) { // this function gets called when simulated
	// first calculate effort
	for (auto &method : *_methods) { method.initializeEffort(); }

	auto detectionProbs = _readDetectionProbFile();
	auto N0             = _readN0FromCommandLine();

	// do we need to simulate covariates?
	bool simulateCovariates = !parameters().parameterExists("simulationFile");

	for (size_t s = 0; s < _epochModel.size(); s++) {
		_epochModel[s]->simulate(*_methods, N0[s], detectionProbs, *_locations, simulateCovariates);
	}
}

//------------------------------------
// Getters

const TUniqueContainer<std::string> &TBirpPrior::locations() const { return *_locations; }

const std::vector<TMethods> &TBirpPrior::methods() const { return *_methods; }

const TUniqueContainer<size_t> &TBirpPrior::timepoints() const { return *_timepoints; }
