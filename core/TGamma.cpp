//
// Created by caduffm on 4/17/23.
//

#include "TGamma.h"

#include "TData.h"
#include "TModels.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TProbabilityDistributions.h"

using namespace coretools::instances;

//------------------------------------------
// TGamma
//------------------------------------------

TGamma::TGamma(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1> *TimesOfChange)
    : _gamma(Gamma), _timesOfChange(TimesOfChange) {
	// find out how many epochs there are
	if (parameters().parameterExists("timesOfChange")) {
		std::vector<std::string> tmp;
		parameters().fillParameterIntoContainer("timesOfChange", tmp, ',');
		_numEpochs = tmp.size() + 1;
	} else {
		_numEpochs = parameters().getParameterWithDefault("numEpochs", 1);
	}
	// default: do not infer times of change
	if (!parameters().parameterExists("inferTimesOfChange")) { _timesOfChange->getDefinition().update(false); }
}

void TGamma::initialize(const TUniqueContainer<size_t> &Timepoints) {
	_timesOfChange->initStorage({_numEpochs - 1}, {std::make_shared<coretools::TNamesIndices>()});
	_gamma->initStorage({_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});

	_fillRho(Timepoints);
}

void TGamma::estimateInitialGamma(const TData &Data, const TModelBase &Model) {
	// TODO: something more sophisticated?
	for (size_t g = 0; g < _gamma->totalSize(); ++g) { _gamma->setInitial(g, 0.0); }

	// set _oldLL and _oldJP
	_curJP = Model.calculateJeffreyPrior(Data, *this);
	_tryJP = _curJP;
}

std::vector<size_t> TGamma::getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const {
	// return start of each epoch, including for a hypothetical epoch after the last epoch
	assert(!Timepoints.empty());

	std::vector<size_t> epochStarts = {0};
	if (_timesOfChange->totalSize() != 0) {
		size_t maxTimeOfChange = _timesOfChange->max_value();
		size_t maxTime         = std::max(Timepoints.back(), maxTimeOfChange + 1);
		for (size_t i = 0; i < _timesOfChange->totalSize(); i++) { epochStarts.push_back(_timesOfChange->value(i)); }
		epochStarts.push_back(maxTime);
	} else {
		epochStarts.push_back(*max_element(Timepoints.cbegin(), Timepoints.cend()));
	}
	return epochStarts;
}

void TGamma::_fillRho(const TUniqueContainer<size_t> &Timepoints) {
	std::vector<size_t> epochStarts = getEpochStarts(Timepoints);
	_rho.zeros((int)Timepoints.size(), (int)_numEpochs);

	for (size_t t = 0; t < Timepoints.size(); ++t) {
		for (size_t e = 1; e < epochStarts.size(); ++e) {
			if (Timepoints[t] <= epochStarts[e - 1]) {
				_rho(t, e - 1) = 0;
			} else if (epochStarts[e - 1] < Timepoints[t] && Timepoints[t] < epochStarts[e]) {
				_rho(t, e - 1) = Timepoints[t] - epochStarts[e - 1];
			} else {
				_rho(t, e - 1) = epochStarts[e] - epochStarts[e - 1];
			}
		}
	}
}

void TGamma::setSimulatedGammaToZero() {
	for (size_t g = 0; g < _gamma->totalSize(); ++g) { _gamma->setInitial(g, 0.0); }
}

// functions to update gamma

bool TGamma::updateSpecificIndex(size_t m) { return _gamma->updateSpecificIndex(m); }

double TGamma::getLogPriorRatio(const TData &Data, const TModelBase &Model) {
	// -> prior on gamma: Jeffrey prior
	_tryJP = Model.calculateJeffreyPrior(Data, *this);
	return _tryJP - _curJP;
}

bool TGamma::acceptOrReject(double logH, size_t m) {
	if (_gamma->acceptOrReject(logH, m)) {
		_curJP = _tryJP;
		return true;
	}
	return false;
}

// getters
size_t TGamma::numEpochs() const { return _numEpochs; }
size_t TGamma::numTimepoints() const { return _rho.n_rows; }

TypeGamma TGamma::gamma(size_t m) const { return _gamma->value(m); }

size_t TGamma::rho(size_t row, size_t col) const { return _rho(row, col); }

//------------------------------------------
// TStochasticPrior
//------------------------------------------

TStochasticPrior::TStochasticPrior(stattools::TParameterTyped<TypePhi, 2> *Phi,
                                   stattools::TParameterTyped<TypeSigma, 1> *Sigma,
                                   stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                                   stattools::TParameterTyped<TypeTime, 1> *TimesOfChange)
    : _phi(Phi), _sigma(Sigma), _gamma(Gamma, TimesOfChange) {}

void TStochasticPrior::initialize(const TUniqueContainer<size_t> &Timepoints) {
	// initialize gamma and times of change
	_gamma.initialize(Timepoints);
}

void TStochasticPrior::estimateInitialGamma(const TData &Data, const TModelBase &Model) {
	_gamma.estimateInitialGamma(Data, Model);
}

void TStochasticPrior::_fillStochasticRho(const TUniqueContainer<size_t> &Timepoints) {
	std::vector<size_t> epochStarts = _gamma.getEpochStarts(Timepoints);
	_stochasticRho.zeros((int)Timepoints.size(), (int)_gamma.numEpochs());

	for (size_t t = 0; t < Timepoints.size(); ++t) {
		for (size_t e = 1; e < epochStarts.size(); ++e) {
			// todo: implement equation 6 (rho for stochastic model)
		}
	}
}

coretools::StrictlyPositive<double> TStochasticPrior::_calculateDelta(size_t k) const {
	// todo: implement equation 5 (calculate delta form exp(sum(rho * gamma))
	// todo: probably store intermediate delta's
	return 0.0;
}

coretools::StrictlyPositive<double> TStochasticPrior::_calculateDeltaTime(size_t k) const {
	// todo: calculate t_k - t_{k-1} -> pre-calculate store in some matrix
}

coretools::Positive<double> TStochasticPrior::_calculateLogTransitionProbability(size_t j, size_t k) const {
	assert(k > 0);

	// calculate mean of normal distribution
	const auto delta            = _calculateDelta(k);
	const auto phi_j_previous_k = _phi->value(_phi->getIndex({j, k - 1}));
	const auto mean             = phi_j_previous_k * delta;

	// calculate sd of normal distribution
	const auto deltaTime = _calculateDeltaTime(k);
	const auto sd        = _sigma->value() * deltaTime;

	const auto phi_j_k = _phi->value(_phi->getIndex({j, k}));
	coretools::probdist::TNormalDistr distr(mean, sd);
	return distr.logDensity(phi_j_k) + distr.logDensity(-phi_j_k);
}

double TStochasticPrior::calculateLLRatioForUpdateGamma() {
	// likelihood: just go through Markov chain and sum up log transition probabilities
	double LL = 0.0;
	for (size_t j = 0; j < _phi->dimensions()[0]; ++j) {     // loop over locations
		for (size_t k = 1; k < _phi->dimensions()[1]; ++k) { // loop over timepoints
			_tryLogTransitionProbs[j][k] = _calculateLogTransitionProbability(j, k);
			LL += _tryLogTransitionProbs[j][k] - _curLogTransitionProbs[j][k];
		}
	}
	return LL;
}

double TStochasticPrior::calculateLL() const {
	// likelihood: just go through Markov chain and sum up log transition probabilities
	double LL = 0.0;
	for (size_t j = 0; j < _phi->dimensions()[0]; ++j) {     // loop over locations
		for (size_t k = 1; k < _phi->dimensions()[1]; ++k) { // loop over timepoints
			LL += _calculateLogTransitionProbability(j, k);
		}
	}
	return LL;
}

void TStochasticPrior::updateGamma(const TData &Data, const TModelBase &Model) {
	for (size_t m = 0; m < _gamma.numEpochs(); m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma.updateSpecificIndex(m)) {
			const double LLRatio       = calculateLLRatioForUpdateGamma();
			const double logPriorRatio = _gamma.getLogPriorRatio(Data, Model); // Jeffrey prior

			const double logH = LLRatio + logPriorRatio;
			if (_gamma.acceptOrReject(logH, m)) {
				// accepted -> update LL and jeffrey's Prior
				swapTryCur(); // -> _oldLL = newLL;
			}
		}
	}
}

double TStochasticPrior::calculateLogPriorRatio(size_t j, size_t k) {
	_tryLogTransitionProbs[j][k] = _calculateLogTransitionProbability(j, k);
	return _tryLogTransitionProbs[j][k] + _tryLogTransitionProbs[j][k];
}

void TStochasticPrior::swapTryCur() {
	// used when gamma updated (-> all transition probabilities change)
	std::swap(_tryLogTransitionProbs, _curLogTransitionProbs);
}

void TStochasticPrior::swapTryCur(size_t j, size_t k) {
	// used when a phi is updated
	std::swap(_tryLogTransitionProbs[j][k], _curLogTransitionProbs[j][k]);
}

void TStochasticPrior::setSimulatedGammaToZero() { _gamma.setSimulatedGammaToZero(); }
