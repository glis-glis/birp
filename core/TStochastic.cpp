//
// Created by caduffm on 4/13/23.
//

#include "TStochastic.h"

#include "coretools/Main/TParameters.h"
#include "coretools/Math/TProbabilityDistributions.h"
#include <utility>

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

void TGamma::initializeInferred(const TUniqueContainer<size_t> &Timepoints) {
	_timesOfChange->initStorage({_numEpochs - 1}, {std::make_shared<coretools::TNamesIndices>()});
	_gamma->initStorage({_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});

	_fillRho(Timepoints);
}

std::vector<size_t> TGamma::_getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const {
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
	std::vector<size_t> epochStarts = _getEpochStarts(Timepoints);
	_rho.zeros((int)Timepoints.size(), (int)_numEpochs);

	for (size_t t = 0; t < Timepoints.size(); ++t) {
		for (size_t e = 1; e < epochStarts.size(); ++e) {
			// todo: implement equation 6 (rho for stochastic model)
		}
	}
}

coretools::StrictlyPositive<double> TGamma::calculateDelta(size_t k) {
	// todo: implement equation 5 (calculate delta form exp(sum(rho * gamma))
	// todo: probably store intermediate delta's
	return 0.0;
}

coretools::StrictlyPositive<double> TGamma::calculateDeltaTime(size_t k) {
	// todo: calculate t_k - t_{k-1}
}

void TGamma::updateGamma(TPhi &Phi, TLikelihoodCalculator &LikelihoodCalculator) {
	for (size_t m = 0; m < _numEpochs; m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma->updateSpecificIndex(m)) {
			const double LLRatio       = Phi.calculateLLRatioForUpdateGamma();
			const double newJP         = LikelihoodCalculator.calculateJeffreyPrior();
			const double logPriorRatio = newJP - _oldJP; // in log

			const double logH = LLRatio + logPriorRatio;
			if (_gamma->acceptOrReject(logH, m)) {
				// accepted -> update LL and jeffrey's Prior
				Phi.swapTryCur(); // -> _oldLL = newLL;
				_oldJP = newJP;
			}
		}
	}
}

//------------------------------------------
// TPhi
//------------------------------------------

TPhi::TPhi(std::vector<stattools::TParameterTyped<TypePhi, 1> *> Phi, stattools::TParameterTyped<TypeSigma, 1> *Sigma,
           TGamma *Gamma)
    : _phi(std::move(Phi)), _sigma(Sigma), _gamma(Gamma) {}

void TPhi::initializeInferred(const TUniqueContainer<size_t> &Timepoints) {
	// todo: different locations have different number of timepoints
	//  -> store phi tidy (for each location same number of timepoints) or ragged?
	//  -> tidy: would need to know if there is data; but is this ok for HMM? Overhead?
	//  -> ragged: would need to know sturcture of data -> birp-independent?
	// _phi->initStorage({...}, {std::make_shared<coretools::TNamesIndices>()});
	_sigma->initStorage();

	// todo: set first phi (k=0) of each location to one (is never updated)

	// todo: resize _tryLogTransitionProbs and _curLogTransitionProbs to correct size (matches size of phi) and fill
	// them from initial values
}

coretools::Positive<double> TPhi::_calculateLogTransitionProbability(size_t j, size_t k) const {
	assert(k > 0);

	// calculate transition probability of k-1 -> k
	const auto &phi = _phi[j];

	// calculate mean of normal distribution
	const auto delta            = _gamma->calculateDelta(k);
	const auto phi_j_previous_k = phi->value(k - 1);
	const auto mean             = phi_j_previous_k * delta;

	// calculate sd of normal distribution
	const auto deltaTime = _gamma->calculateDeltaTime(k);
	const auto sd        = _sigma->value() * deltaTime;

	const auto phi_j_k = phi->value(k);
	coretools::probdist::TNormalDistr distr(mean, sd);
	return distr.logDensity(phi_j_k) + distr.logDensity(-phi_j_k);
}

double TPhi::calculateLL() const {
	// likelihood: just go through Markov chain and sum up log transition probabilities
	double LL = 0.0;
	for (size_t j = 0; j < _phi.size(); ++j) {
		for (size_t k = 1; k < _phi[j]->totalSize(); ++k) { LL += _calculateLogTransitionProbability(j, k); }
	}
	return LL;
}

double TPhi::calculateLLRatioForUpdateGamma() {
	// likelihood: just go through Markov chain and sum up log transition probabilities
	double LL = 0.0;
	for (size_t j = 0; j < _phi.size(); ++j) {
		for (size_t k = 1; k < _phi[j]->totalSize(); ++k) {
			_tryLogTransitionProbs[j][k] = _calculateLogTransitionProbability(j, k);
			LL += _tryLogTransitionProbs[j][k] - _curLogTransitionProbs[j][k];
		}
	}
	return LL;
}

void TPhi::swapTryCur() {
	// cur = try
	// try = cur
	std::swap(_tryLogTransitionProbs, _curLogTransitionProbs);
}

void TPhi::_updatePhi(size_t j, size_t k, TLikelihoodCalculator &LikelihoodCalculator) {
	const auto &phi = _phi[j];
	// propose a new value
	phi->updateSpecificIndex(k);

	const auto logLeftNeighbour  = _calculateLogTransitionProbability(j, k);
	const auto logRightNeighbour = _calculateLogTransitionProbability(j, k + 1);
	// only need to calculate LL for this location j and timepoint k, but across all methods (if they have data for j
	// and k)
	const auto LLRatio           = LikelihoodCalculator.calculateLLRatio(j, k);
	const double logH            = LLRatio + logLeftNeighbour + logRightNeighbour - _curLogTransitionProbs[j][k] -
	                    _curLogTransitionProbs[j][k + 1];
	if (phi->acceptOrReject(logH, k)) {
		// accepted
		_curLogTransitionProbs[j][k]     = logLeftNeighbour;
		_curLogTransitionProbs[j][k + 1] = logRightNeighbour;
		// todo: tell LikelihoodCalculator that update was accepted??
	}
}

void TPhi::_updatePhi(TLikelihoodCalculator &LikelihoodCalculator) {
	for (size_t j = 0; j < _phi.size(); ++j) { // loop over all locations
		for (size_t k = 1; k < _phi[j]->totalSize(); ++k) {
			// loop over all timepoints (ignore first, since this is always 1)
			_updatePhi(j, k);
		}
	}
}

void TPhi::_updateSigma(TLikelihoodCalculator &LikelihoodCalculator) {
	// propose a new value
	_sigma->update();

	double logTransitionRatio = 0.0;
	double LLRatio            = 0.0;
	for (size_t j = 0; j < _phi.size(); ++j) {
		for (size_t k = 1; k < _phi[j]->totalSize(); ++k) {
			_tryLogTransitionProbs[j][k] = _calculateLogTransitionProbability(j, k);
			logTransitionRatio += _tryLogTransitionProbs[j][k] - _curLogTransitionProbs[j][k];
			LLRatio += LikelihoodCalculator.calculateLLRatio(j, k);
		}
	}

	const double logH = logTransitionRatio + LLRatio;
	if (_sigma->acceptOrReject(logH, 0)) {
		// accepted: update transition probabilities
		swapTryCur();
		// todo: tell LikelihoodCalculator that update was accepted??
	}
}

void TPhi::update(TLikelihoodCalculator &LikelihoodCalculator) {
	_updatePhi(LikelihoodCalculator);
	_updateSigma(LikelihoodCalculator)
}
