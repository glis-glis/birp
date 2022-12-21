
//
// Created by Liam on 23/08/22.
//

#include "TEpochModel.h"
#include "TBirpPrior.h"

using namespace coretools::instances;

TEpochModel::TEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                         stattools::TParameterTyped<TypeTime, 1> *TimesOfChange, const std::string &Filename)
    : _gamma(Gamma), _timesOfChange(TimesOfChange) {}

void TEpochModel::initialize() {
	if (!_timesOfChange->hasFixedInitialValue()) { // if user did not set initial values
		_numEpochs = 1;
	}
	_numEpochs = _timesOfChange->totalSize() + 1;
	_gamma->initStorage({_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});
}

std::vector<size_t> TEpochModel::_getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const {
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

// static void _readTimesOfChange(){}

void TEpochModel::_fillRho(const TUniqueContainer<size_t> &Timepoints) {
	std::vector<size_t> epochStarts = _getEpochStarts(Timepoints);
	_rho.zeros((int)Timepoints.size(), (int)_numEpochs);
	// fill rhos
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

double TEpochModel::_calculateLogPhi_k(size_t k) const {
	double sum = 0.0;
	for (size_t m = 0; m < _numEpochs; m++) { sum += (double)_rho(k, m) * _gamma->value(m); }
	return sum;
}

std::vector<double> TEpochModel::_calculatePhi() const {
	std::vector<double> phi(_rho.n_rows);
	for (size_t k = 0; k < _rho.n_rows; k++) { phi[k] = _calculateLogPhi_k(k); }
	coretools::standardizeZeroMean(phi);
	for (auto &p : phi) { p = exp(p); } // make logPhi to phi
	return phi;
}

arma::mat TEpochModel::_calcFisherInfo(const std::vector<double> &p, const std::vector<double> &sigma, TypeCounts nu) const {
	arma::mat R(p.size(), _numEpochs);
	arma::mat RtP(_numEpochs, p.size());
	for (size_t k = 0; k < p.size(); k++) {
		for (size_t m = 0; m < _numEpochs; m++) {
			double denom       = 0.0;
			double rhoSumSigma = 0.0;
			for (size_t l = 0; l < p.size(); l++) {
				denom += sigma[l];
				rhoSumSigma += ((double)_rho(k, m) - (double)_rho(l, m)) * sigma[l];
			}
			denom     = denom * denom;
			R(k, m)   = sigma[k] * rhoSumSigma / denom;
			RtP(m, k) = R(k, m) * 1.0 / p[k]; // write test for this
		}
	}
	// real actual matrix multiplication
	auto tmp = RtP * R;
	return nu * tmp;
}

std::pair<double, double> TEpochModel::_calculateLLandJP(const TSpecies &Data) const {
	std::vector<double> phi = _calculatePhi();
	double LL               = 0.0;
	arma::mat fisherInfo(_numEpochs,_numEpochs, arma::fill::zeros);
	for (size_t i = 0; i < Data.size(); i++) {
		const auto &lm = Data[i];
		std::vector<double> sigma(lm.size());
		for (size_t k = 0; k < lm.size(); k++) {
			// calculate sigma and fill into p (will normalize later)
			sigma[k] = phi[k] * lm[k].effort;
		}
		// now normalize by row sum!
		std::vector<double> p = sigma;
		coretools::normalize(p);
		// add to LL
		for (size_t k = 0; k < lm.size(); k++) { LL += (double)lm[k].counts * log(p[k]); }
		// add to JP
		fisherInfo += _calcFisherInfo(p, sigma, lm.nu());
	}
	double JP = log(arma::det(fisherInfo));  //sqrt drops out in hastings ratio
	return {LL, JP};
}

void TEpochModel::updateGamma(const TSpecies &Data) {
	for (size_t m = 0; m < _numEpochs; m++) { // loop over epochs and update gamma epoch by epoch
		_gamma->updateSpecificIndex(m);

		const auto [newLL, newJP]   = _calculateLLandJP(Data);
		const double LLRatio = newLL - _oldLL;  // in log
		const double logPriorRatio = newJP - _oldJP; // in log

		const double logH = LLRatio + logPriorRatio;
		if (_gamma->acceptOrReject(m, logH)) {
			// accepted -> update LL and jeffrey's Prior
			_oldLL           = newLL;
			_oldJP = newJP;
		}
	}
}

// TO DO NEXT WEEK:
// Stiflä (Madleina)
// Simulator schribä
// remove data sets with less than two time points with effort > 0





// getters for test functions
size_t TEpochModel::rho(size_t row, size_t col) const { return _rho(row, col); }


