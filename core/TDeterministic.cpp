//
// Created by caduffm on 4/13/23.
//

#include "TDeterministic.h"

using namespace coretools::instances;

//------------------------------------------
// TDeterministic
//------------------------------------------

TDeterministic::TDeterministic(TGamma *Gamma) : _gamma(Gamma) { _phi.resize(_gamma->numTimepoints()); }

double TDeterministic::_calculateLogPhi_k(size_t k) const {
	double sum = 0.0;
	for (size_t m = 0; m < _gamma->numEpochs(); m++) { sum += (double)_gamma->rho(k, m) * _gamma->gamma(m); }
	return sum;
}

void TDeterministic::_fillPhi() {
	for (size_t k = 0; k < _gamma->numTimepoints(); k++) { _phi[k] = _calculateLogPhi_k(k); }
	coretools::standardizeZeroMean(_phi);
	for (auto &p : _phi) { p = exp(p); } // make logPhi to phi
}

TypePhi TDeterministic::phi(size_t /*location_id*/, size_t timepoint_id) const { return _phi[timepoint_id]; }

void TDeterministic::updateGamma(const TData &Data, TLikelihood &Likelihood) {
	for (size_t m = 0; m < _gamma->numEpochs(); m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma->updateSpecificIndex(m)) {
			const auto LLRatio         = Likelihood.calculateLLRatio(Data, *this, *_distr);
			const double logPriorRatio = _gamma->getLogPriorRatio(Data, Likelihood, *this); // Jeffrey prior

			const double logH = LLRatio + logPriorRatio;
			if (_gamma->acceptOrReject(logH, m)) { Likelihood.swapTryCur(); }
		}
	}
}