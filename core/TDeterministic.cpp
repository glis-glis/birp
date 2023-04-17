//
// Created by caduffm on 4/13/23.
//

#include "TDeterministic.h"

using namespace coretools::instances;

//------------------------------------------
// TDeterministic
//------------------------------------------

TDeterministic::TDeterministic(TGamma *Gamma) : _gamma(Gamma) {}

double TDeterministic::_calculateLogPhi_k(size_t k) const {
	double sum = 0.0;
	for (size_t m = 0; m < _gamma->numEpochs(); m++) { sum += (double)_gamma->rho(k, m) * _gamma->gamma(m); }
	return sum;
}

std::vector<double> TDeterministic::calculatePhi(size_t /*j*/) const {
	std::vector<double> phi(_gamma->numTimepoints());
	for (size_t k = 0; k < _gamma->numTimepoints(); k++) { phi[k] = _calculateLogPhi_k(k); }
	coretools::standardizeZeroMean(phi);
	for (auto &p : phi) { p = exp(p); } // make logPhi to phi
	return phi;
}

void TDeterministic::updateGamma(TLikelihood &Likelihood) {
	for (size_t m = 0; m < _gamma->numEpochs(); m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma->updateSpecificIndex(m)) {
			const auto LLRatio = Likelihood.calculateLLRatio();
			const double logPriorRatio = _gamma->getLogPriorRatio(Likelihood); // Jeffrey prior

			const double logH = LLRatio + logPriorRatio;
			if (_gamma->acceptOrReject(logH, m)) {
				// todo: tell Likelihood to update likelihood
			}
		}
	}
}