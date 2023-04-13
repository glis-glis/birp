//
// Created by caduffm on 4/13/23.
//

#include "TDeterministic.h"

using namespace coretools::instances;

//------------------------------------------
// TDeterministic
//------------------------------------------

TDeterministic::TDeterministic(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                               stattools::TParameterTyped<TypeTime, 1> *TimesOfChange)
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

void TDeterministic::initializeInferred(const TUniqueContainer<size_t> &Timepoints) {
	_timesOfChange->initStorage({_numEpochs - 1}, {std::make_shared<coretools::TNamesIndices>()});
	_gamma->initStorage({_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});

	_fillRho(Timepoints);
}

void TDeterministic::estimateInitialPriorParameters(TLikelihoodCalculator &LikelihoodCalculator) {
	// TODO: something more sophisticated?
	for (size_t g = 0; g < _gamma->totalSize(); ++g) { _gamma->set(g, 0.0); }

	// set _oldLL and _oldJP
	_oldJP = LikelihoodCalculator.calculateJeffreyPrior();
}

std::vector<size_t> TDeterministic::_getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const {
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

void TDeterministic::_fillRho(const TUniqueContainer<size_t> &Timepoints) {
	std::vector<size_t> epochStarts = _getEpochStarts(Timepoints);
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

double TDeterministic::_calculateLogPhi_k(size_t k) const {
	double sum = 0.0;
	for (size_t m = 0; m < _numEpochs; m++) { sum += (double)_rho(k, m) * _gamma->value(m); }
	return sum;
}

std::vector<double> TDeterministic::calculatePhi(size_t /*j*/) const {
	std::vector<double> phi(_rho.n_rows);
	for (size_t k = 0; k < _rho.n_rows; k++) { phi[k] = _calculateLogPhi_k(k); }
	coretools::standardizeZeroMean(phi);
	for (auto &p : phi) { p = exp(p); } // make logPhi to phi
	return phi;
}

void TDeterministic::updateGamma(TLikelihood &Likelihood) {
	for (size_t m = 0; m < _numEpochs; m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma->updateSpecificIndex(m)) {
			const auto [LLRatio, newJP] = Likelihood.calculateLLRatioandJP();
			const double logPriorRatio  = newJP - _oldJP; // in log

			const double logH = LLRatio + logPriorRatio;
			if (_gamma->acceptOrReject(logH, m)) {
				// accepted -> update LL and jeffrey's Prior
				// todo: tell Likelihood to update likelihood
				_oldJP = newJP;
			}
		}
	}
}

void TDeterministic::_setSimulatedGammaToZero() {
	for (size_t g = 0; g < _gamma->totalSize(); ++g) { _gamma->setInitial(g, 0.0); }
}

// getters for test functions
size_t TDeterministic::rho(size_t row, size_t col) const { return _rho(row, col); }
