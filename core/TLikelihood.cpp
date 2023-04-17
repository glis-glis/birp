//
// Created by Liam on 23/08/22.
//
#include "TLikelihood.h"
#include "TDistributions.h"
#include "TStochastic.h"

//---------------------------------------
// TLikelihood
//---------------------------------------

TLikelihood::TLikelihood(size_t SpeciesIDinUniqueContainer, const TData &Data)
    : _speciesIDinUniqueContainer(SpeciesIDinUniqueContainer) {

	// LL is per method/location combination
	_curLL.resize(Data.numMethLoc(_speciesIDinUniqueContainer), 0.0);
	_tryLL.resize(Data.numMethLoc(_speciesIDinUniqueContainer), 0.0);

	// todo: fill likelihoods after initialization
}

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
				std::vector sigma     = calculateSigma(method[l], Model);
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
                                     const TDistributionBase &Distribution) {
	// loop over all datapoints (methods - locations - timepoints) of a single species
	// used in update of gamma of deterministic model and sigma in stochastic model
	size_t ix_linNumMethLoc = 0;
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
			for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
				std::vector sigma     = calculateSigma(method[l], Model);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to Likelihood
				_tryLL[ix_linNumMethLoc] =
				    Distribution.calculateLL(sigma, method[l], p, index_in_counts, ix_linNumMethLoc);
			}
		}
	}
	return coretools::containerSum(_tryLL) - coretools::containerSum(_curLL);
};

double TLikelihood::calculateLLRatio_perLocation(size_t j, const TData &Data, const TModelBase &Model,
                                                 const TDistributionBase &Distribution) {
	// loop over all methods and timepoints of a particular location of a single species
	// used in update of phi in stochastic model
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
			size_t index_in_counts  = method.getIndexInCounts(_speciesIDinUniqueContainer);
			std::vector sigma       = calculateSigma(method[j], Model);
			std::vector<double> p   = sigma;
			// now normalize by row sum!
			coretools::normalize(p);
			// add to Likelihood
			_tryLL[ix_linNumMethLoc] = Distribution.calculateLL(sigma, method[j], p, index_in_counts, ix_linNumMethLoc);
		}
	}
	return coretools::containerSum(_tryLL) - coretools::containerSum(_curLL);
};

double TLikelihood::calculateLLRatio_perMethod(size_t i, const TData &Data, const TModelBase &Model,
                                               const TDistributionBase &Distribution) {
	// loop over all locations and timepoints of a particular method of a single species
	const auto &method = Data[i];

	if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
		return coretools::containerSum(_curLL);
	}                                            // chame gd hei ga

	for (size_t j = 0; j < method.size(); ++j) { // loop over locations
		const auto &loc        = method[j];
		size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
		std::vector sigma      = calculateSigma(loc, Model);
		std::vector<double> p  = sigma;
		// now normalize by row sum!
		coretools::normalize(p);
		// calculate Likelihood
		size_t ix_linNumMethLoc  = Data.i_j_to_linear(_speciesIDinUniqueContainer, method.method_id(),
		                                              j); // method_id should always be the same as i!
		_tryLL[ix_linNumMethLoc] = Distribution.calculateLL(sigma, method[j], p, index_in_counts, ix_linNumMethLoc);
	}
	return coretools::containerSum(_tryLL);
}

void TLikelihood::swapTryCur() {
	// used when all datapoints are affected (update of deterministic gamma or sigma)
	std::swap(_tryLL, _curLL);
}

void TLikelihood::swapTryCur_perLocation(size_t j, const TData &Data) {
	// used when a location has changed
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
			std::swap(_tryLL[ix_linNumMethLoc], _curLL[ix_linNumMethLoc]);
		}
	}
}

void TLikelihood::swapTryCur_perMethod(size_t i, const TData &Data) {
	// used when a method has changed
	const auto &method = Data[i];
	if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return; }
	for (size_t j = 0; j < method.size(); ++j) { // loop over locations
		size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, method.method_id(), j);
		std::swap(_tryLL[ix_linNumMethLoc], _curLL[ix_linNumMethLoc]);
	}
}

void TLikelihood::setLLs(const std::vector<double> &LLs) { _oldLL = LLs; }
const std::vector<double> &TLikelihood::oldLLs() const { return _oldLL; }

// TODO
// automatically remove data sets with less than two time points with effort > 0?
