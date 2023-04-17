//
// Created by Liam on 23/08/22.
//
#include "TLikelihood.h"
#include "TDistributions.h"
#include "TStochastic.h"

//---------------------------------------
// TLikelihood
//---------------------------------------

TLikelihood::TLikelihood(size_t SpeciesIDinUniqueContainer) : _speciesIDinUniqueContainer(SpeciesIDinUniqueContainer) {}

// todo: also standardizeZeroMean for stochastic phi?

std::vector<double> TLikelihood::calculateSigma(const TLocations &loc, const TModelBase &Model) {
	std::vector<double> sigma(loc.size());
	for (size_t k = 0; k < loc.size(); k++) {
		// calculate sigma and fill into p (will normalize later)
		sigma[k] = Model.phi(loc.location_id(), loc[k].timepoint_id()) * loc[k].effort();
	}
	return sigma;
}

arma::mat TLikelihood::_calcFisherInfo(const TGamma &Gamma, const std::vector<double> &p,
                                       const std::vector<double> &sigma, TypeNu nu) {
	arma::mat R(p.size(), Gamma.numEpochs());
	arma::mat RtP(Gamma.numEpochs(), p.size());
	for (size_t k = 0; k < p.size(); k++) {
		for (size_t m = 0; m < Gamma.numEpochs(); m++) {
			double denom       = 0.0;
			double rhoSumSigma = 0.0;
			for (size_t l = 0; l < p.size(); l++) {
				denom += sigma[l];
				rhoSumSigma += ((double)Gamma.rho(k, m) - (double)Gamma.rho(l, m)) * sigma[l];
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

double TLikelihood::calculateJeffreyPrior(const TData &Data, const TGamma &Gamma, const TModelBase &Model) const {
	// used in update of gamma of stochastic model
	arma::mat fisherInfo(Gamma.numEpochs(), Gamma.numEpochs(), arma::fill::zeros);
	size_t ix_linNumMethLoc = 0;
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
			for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
				std::vector sigma     = _calcSigma(method[l], Model);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to JP
				fisherInfo += _calcFisherInfo(Gamma, p, sigma, method[l].nu(index_in_counts));
			}
		}
	}
	return log(arma::det(fisherInfo)); // sqrt drops out in hastings ratio
}

double TLikelihood::calculateLLRatio(const TData &Data, const TModelBase &Model,
                                                  const TDistributionBase &Distribution) const {
	// used in update of gamma of deterministic model and sigma in stochastic model
	std::vector<double> LL(Data.numMethLoc(_speciesIDinUniqueContainer), 0.0); // LL is per method/location combination
	size_t ix_linNumMethLoc = 0;
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
			for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
				std::vector sigma     = _calcSigma(method[l], Model);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to Likelihood
				LL[ix_linNumMethLoc] = Distribution.calculateLL(sigma, method[l], p, index_in_counts, ix_linNumMethLoc);
			}
		}
	}
	return coretools::containerSum(LL); // todo: vector or single number (sum)? Also: currently not a ratio
};

double TLikelihood::calculateLLRatio(size_t j, const TData &Data, const TModelBase &Model,
                                     const TDistributionBase &Distribution) const {
	// used in update of phi in stochastic model
	double LL = 0.0;
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
			size_t index_in_counts  = method.getIndexInCounts(_speciesIDinUniqueContainer);
			std::vector sigma       = _calcSigma(method[j], Model);
			std::vector<double> p   = sigma;
			// now normalize by row sum!
			coretools::normalize(p);
			// add to Likelihood
			LL += Distribution.calculateLL(sigma, method[j], p, index_in_counts, ix_linNumMethLoc);
		}
	}
	return LL; // todo: vector or single number (sum)? Also: currently not a ratio
};

std::vector<double> TLikelihood::calculateLLForBetaUpdate(size_t i, const TData &Data, const TModelBase &Model,
                                                          const TDistributionBase &Distribution) const {
	std::vector<double> LL = _oldLL; // LL is per method/location combination
	const auto &method     = Data[i];

	if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return LL; } // chame gd hei ga

	for (size_t j = 0; j < method.size(); ++j) {                                 // loop over locations
		const auto &loc        = method[j];
		size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
		std::vector sigma      = _calcSigma(loc, Model);
		std::vector<double> p  = sigma;
		// now normalize by row sum!
		coretools::normalize(p);
		// calculate Likelihood
		size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, method.method_id(), j); // method_id should always be the same as i!
		LL[ix_linNumMethLoc]    =Distribution.calculateLL(sigma, method[j], p, index_in_counts, ix_linNumMethLoc);
	}
	return LL;
}

void TLikelihood::setLLs(const std::vector<double> &LLs) { _oldLL = LLs; }
const std::vector<double> &TLikelihood::oldLLs() const { return _oldLL; }

// TODO
// automatically remove data sets with less than two time points with effort > 0?
