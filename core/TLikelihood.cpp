//
// Created by Liam on 23/08/22.
//
#include "TLikelihood.h"
#include "TDistributions.h"
#include "TStochastic.h"

//---------------------------------------
// TLikelihood
//---------------------------------------

TLikelihood::TLikelihood() {}

/*
// todo: also standardizeZeroMean for stochastic phi?
double TLikelihood::_calculateLogPhi_k(size_t k) const {
    double sum = 0.0;
    for (size_t m = 0; m < _numEpochs; m++) { sum += (double)_rho(k, m) * _gamma->value(m); }
    return sum;
}

std::vector<double> TLikelihood::_calculatePhi() const {
    std::vector<double> phi(_rho.n_rows);
    for (size_t k = 0; k < _rho.n_rows; k++) { phi[k] = _calculateLogPhi_k(k); }
    coretools::standardizeZeroMean(phi);
    for (auto &p : phi) { p = exp(p); } // make logPhi to phi
    return phi;
}
*/

arma::mat TLikelihood::_calcFisherInfo(const std::vector<double> &p, const std::vector<double> &sigma,
                                       TypeNu nu) const {
	arma::mat R(p.size(), _numEpochs);
	arma::mat RtP(_numEpochs, p.size());
	for (size_t k = 0; k < p.size(); k++) {
		for (size_t m = 0; m < _numEpochs; m++) {
			double denom       = 0.0;
			double rhoSumSigma = 0.0;
			for (size_t l = 0; l < p.size(); l++) {
				denom += sigma[l];
				rhoSumSigma += ((double)_rho(k, m) - (double)_rho(l, m)) *
				               sigma[l]; // todo: needs access to "standard" rho -> how to get this?
			}
			denom     = denom * denom;
			R(k, m)   = sigma[k] * rhoSumSigma / denom;
			RtP(m, k) = R(k, m) * 1.0 / p[k];
		}
	}
	// real actual matrix multiplication
	auto tmp = RtP * R;
	return nu * tmp;
}

std::vector<double> TLikelihood::_calcSigma(const TLocations &loc, const std::vector<double> &phi) const {
	std::vector<double> sigma(loc.size());
	for (size_t k = 0; k < loc.size();
	     k++) { // todo: is this ok? one location might not have all timepoints -> need to take phi of correct timepoint
		// calculate sigma and fill into p (will normalize later)
		sigma[k] = phi[k] * loc[k].effort();
	}
	return sigma;
}

double TLikelihood::calculateJeffreyPrior(const TData &Data, size_t SpeciesIDinUniqueContainer,
                                          const TStochastic &Stochastic) const {
	// used in update of gamma of stochastic model
	arma::mat fisherInfo(_numEpochs, _numEpochs, arma::fill::zeros);
	size_t ix_linNumMethLoc = 0;
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(SpeciesIDinUniqueContainer)) {
			size_t index_in_counts = method.getIndexInCounts(SpeciesIDinUniqueContainer);
			for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
				std::vector<double> phi =
				    Stochastic.calculatePhi(l); // todo: is l the correct index? Not every method has every location!
				std::vector sigma     = _calcSigma(method[l], phi);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to JP
				fisherInfo += _calcFisherInfo(p, sigma, method[l].nu(index_in_counts));
			}
		}
	}
	return log(arma::det(fisherInfo)); // sqrt drops out in hastings ratio
}

double TLikelihood::calculateLLRatio(const TData &Data, size_t j, size_t k, size_t SpeciesIDinUniqueContainer,
                                     TDistributionBase *Distribution, const TStochastic &Stochastic) const {
	// todo: how to have less code duplication? (first lines are the same for all LL functions)
	// used in update of phi of stochastic model
	double LL               = 0.0;
	size_t ix_linNumMethLoc = 0;
	for (size_t i = 0; i < Data.size(); ++i) { // loop over all methods
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(SpeciesIDinUniqueContainer)) {
			size_t index_in_counts = method.getIndexInCounts(SpeciesIDinUniqueContainer);
			for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
				if (l != j) { continue; }
				std::vector<double> phi =
				    Stochastic.calculatePhi(l); // todo: is l the correct index? Not every method has every location!
				std::vector sigma     = _calcSigma(method[l], phi);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to Likelihood
				LL += Distribution->calculateLL(sigma, method[l], p, index_in_counts, ix_linNumMethLoc);
			}
		}
	}
	return LL;
}

double TLikelihood::calculateLL(const std::vector<TMethods> &Methods) const {
	const auto [newLL, newJP] = _calculateLLandJP(Methods);
	return coretools::containerSum(newLL);
}

std::vector<double> TLikelihood::calculateLLForBetaUpdate(const TMethods &Method) const {
	std::vector<double> LL = _oldLL; // LL is per method/location combination

	if (!Method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return LL; } // chame gd hei ga

	std::vector<double> phi = _calculatePhi();
	for (size_t j = 0; j < Method.size(); ++j) { // loop over locations
		const auto &loc        = Method[j];
		size_t index_in_counts = Method.getIndexInCounts(_speciesIDinUniqueContainer);
		std::vector sigma      = _calcSigma(loc, phi);
		std::vector<double> p  = sigma;
		// now normalize by row sum!
		coretools::normalize(p);
		// calculate Likelihood
		size_t ij = _i_j_to_linear[Method.method_id()][j];
		LL[ij]    = _calculateLL(sigma, loc, p, index_in_counts, ij);
	}
	return LL;
}

void TLikelihood::setLLs(const std::vector<double> &LLs) { _oldLL = LLs; }
const std::vector<double> &TLikelihood::oldLLs() const { return _oldLL; }

// TODO
// automatically remove data sets with less than two time points with effort > 0?
