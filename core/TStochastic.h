//
// Created by caduffm on 4/13/23.
//

#ifndef BIRP_TSTOCHASTIC_H
#define BIRP_TSTOCHASTIC_H

#include "BirpTypes.h"
#include "TData.h"
#include "TDistributions.h"
#include "TLikelihood.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>
#include "TGamma.h"

//------------------------------------------
// TStochasticPrior
//------------------------------------------

class TStochasticPrior {
	// yellow box
protected:
	stattools::TParameterTyped<TypePhi, 2> *_phi     = nullptr;
	stattools::TParameterTyped<TypeSigma, 1> *_sigma = nullptr;

	// non-owning pointer to TGamma
	TGamma *_gamma = nullptr;

	// rho
	coretools::TMatrix<size_t> _stochasticRho;

	// store log transition probabilities: per location and timepoint
	std::vector<std::vector<coretools::Positive<double>>> _tryLogTransitionProbs;
	std::vector<std::vector<coretools::Positive<double>>> _curLogTransitionProbs;

	// fill rho
	void _fillStochasticRho(const TUniqueContainer<size_t> &Timepoints);
	coretools::StrictlyPositive<double> _calculateDelta(size_t k) const;
	coretools::StrictlyPositive<double> _calculateDeltaTime(size_t k) const;

	// transition probabilities
	coretools::Positive<double> _calculateLogTransitionProbability(size_t j, size_t k) const;

public:
	TStochasticPrior(stattools::TParameterTyped<TypePhi, 2> * Phi,
	                 stattools::TParameterTyped<TypeSigma, 1> *Sigma, TGamma *Gamma);
	~TStochasticPrior() = default;

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints);

	// update gamma
	void updateGamma(TLikelihood &Likelihood);
	double calculateLL() const;
	double calculateLLRatioForUpdateGamma();
	void swapTryCur();
	void swapTryCur(size_t j, size_t k);

	double calculateLogPriorRatio(size_t j, size_t k);
};

//------------------------------------------
// TStochastic
//------------------------------------------

class TStochastic {
	// green box
protected:
	// parameters
	stattools::TParameterTyped<TypePhi, 2> * _phi;
	stattools::TParameterTyped<TypeSigma, 1> *_sigma = nullptr;

	// ptr to prior above
	TStochasticPrior *_prior;

	TDistributionBase *_distr; // Poisson or negative binomial

	// update functions
	void _updatePhi(size_t j, size_t k, TLikelihood &Likelihood);
	void _updateSigma(TLikelihood &Likelihood);
	void _updatePhi(TLikelihood &Likelihood);

public:
	TStochastic(stattools::TParameterTyped<TypePhi, 2> *Phi,
	            stattools::TParameterTyped<TypeSigma, 1> *Sigma, TStochasticPrior *Prior);
	virtual ~TStochastic() = default;

	void update(TLikelihood &Likelihood);

	coretools::StrictlyPositive<double> phi(size_t j, size_t k) const;
};

#endif // BIRP_TSTOCHASTIC_H
