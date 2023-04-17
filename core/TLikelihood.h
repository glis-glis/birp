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

class TGamma;
class TDistributionBase;

//---------------------------------------------------
// TModelBase
// pure abstract base class of TStochastic and TDeterministic
//---------------------------------------------------

class TModelBase {
public:
	virtual ~TModelBase()                                             = default;
	// abstract function for phi
	virtual TypePhi phi(size_t location_id, size_t timepoint_k) const = 0;
};

//---------------------------------------------------
// TLikelihood
//---------------------------------------------------

class TLikelihood {
	// the flower
private:
	size_t _speciesIDinUniqueContainer = 0;

	// temporary values, stored for speed
	std::vector<double> _oldLL;

	// likelihood
	static std::vector<double> _calcSigma(const TLocations &loc, const TModelBase &Model);
	static arma::mat _calcFisherInfo(const TGamma &Gamma, const std::vector<double> &p,
	                                 const std::vector<double> &sigma, TypeNu nu);

protected:
public:
	TLikelihood(size_t SpeciesIDinUniqueContainer);
	virtual ~TLikelihood() = default;

	// calculate sigma
	static std::vector<double> calculateSigma(const TLocations &loc, const TModelBase &Model);

	double calculateJeffreyPrior(const TData &Data, const TGamma &Gamma, const TModelBase &Model) const;
	double calculateLLRatio(const TData &Data, const TModelBase &Model, const TDistributionBase &Distribution) const;
	double calculateLLRatio(size_t j, const TData &Data, const TModelBase &Model,
	                        const TDistributionBase &Distribution) const;

	std::vector<double> calculateLLForBetaUpdate(size_t i, const TData &Data, const TModelBase &Model,
	                                             const TDistributionBase &Distribution) const;
	// get and set likelihoods
	void setLLs(const std::vector<double> &LLs);
	const std::vector<double> &oldLLs() const;
};

#endif // BIRP_TLIKELIHOOD_H
