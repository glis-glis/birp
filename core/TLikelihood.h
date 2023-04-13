//
// Created by caduffm on 4/13/23.
//

#ifndef BIRP_TLIKELIHOOD_H
#define BIRP_TLIKELIHOOD_H

#include "BirpTypes.h"
#include "TData.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>

class TStochastic; // forward declaration

//---------------------------------------
// TLikelihood
// for both stochastic and deterministic!
//---------------------------------------

class TLikelihood {
	// the flower
protected:
	// temporary values, stored for speed
	std::vector<double> _oldLL;

	// likelihood
	std::vector<double> _calcSigma(const TLocations &loc, const std::vector<double> &phi) const;
	arma::mat _calcFisherInfo(const std::vector<double> &p, const std::vector<double> &sigma, TypeNu nu) const;

public:
	TLikelihood();
	virtual ~TLikelihood() = default;

	virtual void estimateInitialPriorParameters(const TData &Methods);

	template<typename TypePhiCalculator>
	std::pair<double, double> calculateLLRatioandJP(const TData &Data, size_t SpeciesIDinUniqueContainer,
	                                                TDistributionBase *Distribution,
	                                                const TypePhiCalculator &PhiCalculator) const {
		// used in update of gamma of deterministic model
		std::vector<double> LL(Data.numMethLoc(SpeciesIDinUniqueContainer),
		                       0.0); // LL is per method/location combination
		arma::mat fisherInfo(_numEpochs, _numEpochs, arma::fill::zeros);
		size_t ix_linNumMethLoc = 0;
		for (size_t i = 0; i < Data.size(); ++i) {
			const auto &method = Data[i];
			if (method.hasDataForSpeciesID(SpeciesIDinUniqueContainer)) {
				size_t index_in_counts = method.getIndexInCounts(SpeciesIDinUniqueContainer);
				for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
					const std::vector<double> &phi = PhiCalculator.calculatePhi(
					    l); // todo: is l the correct index? Not every method has every location!
					std::vector sigma     = _calcSigma(method[l], phi);
					std::vector<double> p = sigma;
					// now normalize by row sum!
					coretools::normalize(p);
					// add to Likelihood
					LL[ix_linNumMethLoc] =
					    Distribution->calculateLL(sigma, method[l], p, index_in_counts, ix_linNumMethLoc);
					// add to JP
					fisherInfo += _calcFisherInfo(p, sigma, method[l].nu(index_in_counts));
				}
			}
		}
		double JP = log(arma::det(fisherInfo)); // sqrt drops out in hastings ratio
		return {LL, JP};
	};

	double calculateJeffreyPrior(const TData &Data, size_t SpeciesIDinUniqueContainer,
	                             const TStochastic &Stochastic) const;
	double calculateLLRatio(const TData &Data, size_t j, size_t k, size_t SpeciesIDinUniqueContainer,
	                        TDistributionBase *Distribution, const TStochastic &Stochastic) const;

	std::vector<double> calculateLLForBetaUpdate(const TMethods &Method) const;

	// get and set likelihoods
	void setLLs(const std::vector<double> &LLs);
	const std::vector<double> &oldLLs() const;
};

#endif // BIRP_TLIKELIHOOD_H
