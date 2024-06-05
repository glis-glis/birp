
#include "TBirpPrior.h"
#include "TBirpCore.h"
#include "TModels.h"
#include "coretools/Storage/TDataFile.h"
#include <utility>

using namespace coretools::instances;

TBirpPrior::TBirpPrior(TData &Data, const TUniqueContainer<std::string> &Locations,
                       const TUniqueContainer<TypeTime> &TimePoints,
                       const TUniqueContainer<std::string> &CovariateEffortNames,
                       const TUniqueContainer<std::string> &CovariateDetectionNames, size_t numSpecies,
                       TBirpCommonModel &CommonModel, std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn,
                       std::unique_ptr<TBirpStochasticAddOn> &StochasticAddOn, std::string Prefix)
    : _data(Data), _locations(Locations), _timepoints(TimePoints), _covariateEffortNames(CovariateEffortNames),
      _covariateDetectionNames(CovariateDetectionNames), _numSpecies(numSpecies), _out(std::move(Prefix)) {

	// stattools: tell which ones are prior parameters
	addPriorParameter({&CommonModel.gamma, &CommonModel.timesOfChange});
	for (auto &a : CommonModel.alpha) { addPriorParameter(a.get()); }
	for (auto &b : CommonModel.beta0) { addPriorParameter(b.get()); }
	for (auto &b : CommonModel.beta) { addPriorParameter(b.get()); }
	if (_numSpecies > 1) { UERROR("Birp only supports single species at the moment.."); }
	if (NegBinomAddOn) { addPriorParameter({&NegBinomAddOn->muOrN, &NegBinomAddOn->b}); }
	if (StochasticAddOn) { addPriorParameter({&StochasticAddOn->logPhi, &StochasticAddOn->logSigma}); }

	// create models
	size_t speciesID = 0;
	if (StochasticAddOn) {
		_models.emplace_back(std::make_unique<TStochastic>(
		    speciesID, Data, Locations.size(), TimePoints.size(), NegBinomAddOn, &CommonModel.gamma,
		    &CommonModel.timesOfChange, &StochasticAddOn->logPhi, &StochasticAddOn->logSigma));
	} else {
		_models.emplace_back(std::make_unique<TDeterministic>(speciesID, Data, NegBinomAddOn, &CommonModel.gamma,
		                                                      &CommonModel.timesOfChange));
	}

	for (size_t i = 0; i < Data.size(); ++i) {
		_data[i].addAlphaBeta(CommonModel.alpha[i], CommonModel.beta0[i], CommonModel.beta[i]);
	}
}

std::string TBirpPrior::name() const { return "birp"; }

//------------------------------------
// Initialization

void TBirpPrior::initialize() {
	// set size of parameters etc.
	for (auto &_method : _data) { _method.initialize(_covariateEffortNames, _covariateDetectionNames, this); }
	for (size_t s = 0; s < _numSpecies; s++) { _models[s]->initialize(_data, _locations, _timepoints, this); }
}

void TBirpPrior::guessInitialValues() {
	// set initial values of parameters
	for (auto &_method : _data) { _method.estimateInitialAlphaBeta(); }
	for (size_t s = 0; s < _numSpecies; ++s) { _models[s]->estimateInitial(_data, _timepoints); }

	// now initialize temporary variables (LL, JP etc)
	for (size_t s = 0; s < _numSpecies; ++s) { _models[s]->initializeTempVariables(_data); }
}

//------------------------------------
// MCMC

void TBirpPrior::updateParams() {
	for (size_t s = 0; s < _numSpecies; ++s) { _models[s]->update(_data); }
	for (auto &_method : _data) { _method.updateAlphaBeta(_data, _models); }
}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {
	// calculates likelihood
	double sum = 0.0;
	for (size_t s = 0; s < _numSpecies; ++s) { sum += _models[s]->sumCurLL(); }
	return sum;
}

void TBirpPrior::burninHasFinished() {
	for (auto &model : _models) { model->resetSummaryGammaPosterior(); }
}

void TBirpPrior::writeSummaryGammaPosterior(const std::string &Prefix) {
	// todo: make sure for multispecies to add name of species to Prefix!
	for (auto &model : _models) { model->writeSummaryGammaPosterior(Prefix); }
}

//------------------------------------
// Simulations

std::vector<double> TBirpPrior::_readNBarFromCommandLine() const {
	std::vector<double> n_bar;
	std::vector<double> n_bar_default = {1000.0};
	parameters().fill("n_bar", n_bar, n_bar_default);
	if (n_bar.size() == 1) {
		n_bar.resize(_data.size(), n_bar[0]); // same for all methods
	}
	if (n_bar.size() != _data.size()) {
		UERROR("Size of n_bar (", n_bar.size(), ") does not match the number of methods (", _data.size(), ")!");
	}
	logfile().list("N_bar = ", n_bar, " (argument 'n_bar').");
	return n_bar;
}

void TBirpPrior::_simulateUnderPrior(Storage *) { // this function gets called when simulated
	// first calculate effort
	for (auto &method : _data) { method.initializeEffort(); }

	auto n_bar = _readNBarFromCommandLine();

	// do we need to simulate covariates?
	bool simulateCovariates = !parameters().exists("simulationFile");

	for (const auto &model : _models) { model->simulate(_data, n_bar, simulateCovariates); }
}

//------------------------------------
// Getters

const TUniqueContainer<std::string> &TBirpPrior::locations() const { return _locations; }

const TData &TBirpPrior::data() const { return _data; }

const TUniqueContainer<TypeTime> &TBirpPrior::timepoints() const { return _timepoints; }
