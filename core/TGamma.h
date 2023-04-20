//
// Created by caduffm on 4/17/23.
//

#ifndef BIRP_TGAMMA_H
#define BIRP_TGAMMA_H

#include "BirpTypes.h"
#include "TUniqueContainer.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>

class TModelBase;
class TData;

//------------------------------------------
// TGamma
//------------------------------------------

class TGamma {
	// wraps parameters gamma and timesOfChange
protected:
	// stattools parameters
	stattools::TParameterTyped<TypeGamma, 1> *_gamma        = nullptr;
	stattools::TParameterTyped<TypeTime, 1> *_timesOfChange = nullptr;

	// number of epochs
	size_t _numEpochs = 0;

	// rho
	coretools::TMatrix<size_t> _rho;

	// Jeffrey's prior
	double _tryJP = 0.0;
	double _curJP = 0.0;

	void _fillRho(const TUniqueContainer<size_t> &Timepoints);

public:
	TGamma(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	~TGamma() = default;

	void initialize(const TUniqueContainer<size_t> &Timepoints);
	void estimateInitialGamma(const TData &Data, const TModelBase &Model);

	std::vector<size_t> getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const;

	double getLogPriorRatio(const TData &Data, const TModelBase &Model);

	bool updateSpecificIndex(size_t m);
	bool acceptOrReject(double logH, size_t m);

	void setSimulatedGammaToZero();

	size_t numEpochs() const;
	size_t numTimepoints() const;
	TypeGamma gamma(size_t m) const;
	size_t rho(size_t row, size_t col) const;
};

//------------------------------------------
// TStochasticPrior
//------------------------------------------

class TStochasticPrior {
	// yellow box
protected:
	stattools::TParameterTyped<TypePhi, 2> *_phi     = nullptr;
	stattools::TParameterTyped<TypeSigma, 1> *_sigma = nullptr;

	// TGamma
	TGamma _gamma;

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
	TStochasticPrior(stattools::TParameterTyped<TypePhi, 2> *Phi, stattools::TParameterTyped<TypeSigma, 1> *Sigma,
	                 stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	                 stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	~TStochasticPrior() = default;

	void initialize(const TUniqueContainer<size_t> &Timepoints);
	void estimateInitialGamma(const TData &Data, const TModelBase &Model);

	// update gamma
	void updateGamma(const TData &Data, const TModelBase &Model);
	double calculateLL() const;
	double calculateLLRatioForUpdateGamma();
	void swapTryCur();
	void swapTryCur(size_t j, size_t k);

	double calculateLogPriorRatio(size_t j, size_t k);

	void setSimulatedGammaToZero();
};

#endif // BIRP_TGAMMA_H
