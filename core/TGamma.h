//
// Created by caduffm on 4/17/23.
//

#ifndef BIRP_TGAMMA_H
#define BIRP_TGAMMA_H

#include "BirpTypes.h"
#include "TData.h"
#include "TDistributions.h"
#include "TLikelihood.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>

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

	void _setSimulatedGammaToZero();

public:
	TGamma(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	~TGamma() = default;

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints);
	void estimateInitialGamma(TLikelihood &Likelihood, const TData &Data, const TModelBase &Model);

	std::vector<size_t> getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const;

	double getLogPriorRatio(const TData &Data, const TLikelihood &Likelihood, const TModelBase &Model);

	bool updateSpecificIndex(size_t m);
	bool acceptOrReject(double logH, size_t m);

	size_t numEpochs() const;
	size_t numTimepoints() const;
	TypeGamma gamma(size_t m) const;
	size_t rho(size_t row, size_t col) const;
};

#endif // BIRP_TGAMMA_H
