//
// Created by caduffm on 4/13/23.
//

#ifndef BIRP_TSTOCHASTIC_H
#define BIRP_TSTOCHASTIC_H

#include "BirpTypes.h"
#include "TData.h"
#include "TDistributions.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>
#include "TLikelihood.h"

//------------------------------------------
// TGamma
//------------------------------------------

class TPhi; // forward declaration

class TGamma {
	// yellow box
protected:
	// Birp parameters and co-variates
	stattools::TParameterTyped<TypeGamma, 1> *_gamma        = nullptr;
	stattools::TParameterTyped<TypeTime, 1> *_timesOfChange = nullptr;

	// number of epochs
	size_t _numEpochs = 0;

	// rho
	coretools::TMatrix<size_t> _rho;

	// Jeffrey's prior
	double _oldJP = 0.0;

	void _fillRho(const TUniqueContainer<size_t> &Timepoints);
	std::vector<size_t> _getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const;

public:
	TGamma(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	~TGamma() = default;

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints);

	coretools::StrictlyPositive<double> calculateDelta(size_t k);
	coretools::StrictlyPositive<double> calculateDeltaTime(size_t k);

	void updateGamma(TPhi &Phi, TLikelihood &Likelihood);
};

//------------------------------------------
// TPhi
//------------------------------------------

class TPhi {
	// link between yellow and green box
protected:
	std::vector<stattools::TParameterTyped<TypePhi, 1> *> _phi; // vector over all locations
	stattools::TParameterTyped<TypeSigma, 1> *_sigma = nullptr;

	// non-owning pointer to TGamma
	TGamma *_gamma = nullptr;

	// store log transition probabilities: per location and timepoint
	std::vector<std::vector<coretools::Positive<double>>> _tryLogTransitionProbs;
	std::vector<std::vector<coretools::Positive<double>>> _curLogTransitionProbs;

	// transition probabilities
	coretools::Positive<double> _calculateLogTransitionProbability(size_t j, size_t k) const;

	// update functions
	void _updatePhi(size_t j, size_t k, TLikelihood &Likelihood);
	void _updateSigma(TLikelihood &Likelihood);

public:
	TPhi(std::vector<stattools::TParameterTyped<TypePhi, 1> *> Phi, stattools::TParameterTyped<TypeSigma, 1> *Sigma,
	     TGamma *Gamma);
	~TPhi() = default;

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints);

	double calculateLL() const;
	double calculateLLRatioForUpdateGamma();
	void swapTryCur();

	void update(TLikelihood &Likelihood);
};

//------------------------------------------
// TStochastic
//------------------------------------------

class TStochastic {
	// green box
protected:
	TPhi _phi;
	TDistributionBase *_distr; // Poisson or negative binomial

public:
	TStochastic();
	virtual ~TStochastic() = default;
};

#endif // BIRP_TSTOCHASTIC_H
