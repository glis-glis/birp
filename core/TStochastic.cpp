//
// Created by caduffm on 4/13/23.
//

#include "TStochastic.h"

#include "coretools/Main/TParameters.h"
#include "coretools/Math/TProbabilityDistributions.h"
#include <utility>

using namespace coretools::instances;

//------------------------------------------
// TStochasticPrior
//------------------------------------------

TStochasticPrior::TStochasticPrior(stattools::TParameterTyped<TypePhi, 2> *Phi,
                                   stattools::TParameterTyped<TypeSigma, 1> *Sigma, TGamma *Gamma)
    : _phi(Phi), _sigma(Sigma), _gamma(Gamma) {}

void TStochasticPrior::initializeInferred(const TUniqueContainer<size_t> &Timepoints) {
	_gamma->initializeInferred(Timepoints);
}

void TStochasticPrior::_fillStochasticRho(const TUniqueContainer<size_t> &Timepoints) {
	std::vector<size_t> epochStarts = _gamma->getEpochStarts(Timepoints);
	_stochasticRho.zeros((int)Timepoints.size(), (int)_gamma->numEpochs());

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

void TStochasticPrior::updateGamma(const TData &Data, const TLikelihood &Likelihood, const TModelBase &Model) {
	for (size_t m = 0; m < _gamma->numEpochs(); m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma->updateSpecificIndex(m)) {
			const double LLRatio       = calculateLLRatioForUpdateGamma();
			const double logPriorRatio = _gamma->getLogPriorRatio(Data, Likelihood, Model); // Jeffrey prior

			const double logH = LLRatio + logPriorRatio;
			if (_gamma->acceptOrReject(logH, m)) {
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

//------------------------------------------
// TStochastic
//------------------------------------------

TStochastic::TStochastic(stattools::TParameterTyped<TypePhi, 2> *Phi, stattools::TParameterTyped<TypeSigma, 1> *Sigma,
                         TStochasticPrior *Prior)
    : _phi(Phi), _sigma(Sigma), _prior(Prior) {}

void TStochastic::_updatePhi(size_t j, size_t k, const TData &Data, TLikelihood &Likelihood) {
	// propose a new value
	const size_t ix = _phi->getIndex({j, k});
	_phi->updateSpecificIndex(ix);

	const auto logPriorRatio = _prior->calculateLogPriorRatio(j, k) + _prior->calculateLogPriorRatio(j, k + 1);
	// only need to calculate LL for this location j and timepoint k, but across all methods (if they have data for j
	// and k)
	const auto LLRatio       = Likelihood.calculateLLRatio_perLocation(j, Data, *this, *_distr);
	const double logH        = LLRatio + logPriorRatio;
	if (_phi->acceptOrReject(logH, ix)) {
		_prior->swapTryCur(j, k);
		_prior->swapTryCur(j, k + 1);
		Likelihood.swapTryCur_perLocation(j, Data);
	}
}

void TStochastic::_updatePhi(const TData &Data, TLikelihood &Likelihood) {
	for (size_t j = 0; j < _phi->dimensions()[0]; ++j) { // loop over all locations
		for (size_t k = 1; k < _phi[j].dimensions()[1]; ++k) {
			// loop over all timepoints (ignore first, since this is always 1)
			_updatePhi(j, k, Data, Likelihood);
		}
	}
}

void TStochastic::_updateSigma(const TData &Data, TLikelihood &Likelihood) {
	// propose a new value
	_sigma->update();

	double logTransitionRatio = 0.0;
	double LLRatio            = 0.0;
	for (size_t j = 0; j < _phi->dimensions()[0]; ++j) { // loop over all locations
		for (size_t k = 1; k < _phi[j].dimensions()[1]; ++k) {
			logTransitionRatio += _prior->calculateLogPriorRatio(j, k);
		}
	}
	LLRatio = Likelihood.calculateLLRatio(Data, *this, *_distr);

	const double logH = logTransitionRatio + LLRatio;
	if (_sigma->acceptOrReject(logH, 0)) {
		// accepted: update transition probabilities
		_prior->swapTryCur();
		Likelihood.swapTryCur();
	}
}

void TStochastic::update(const TData &Data, TLikelihood &Likelihood) {
	_updatePhi(Data, Likelihood);
	_updateSigma(Data, Likelihood);
}

coretools::StrictlyPositive<double> TStochastic::phi(size_t j, size_t k) const {
	return _phi->value(_phi->getIndex({j, k}));
}
