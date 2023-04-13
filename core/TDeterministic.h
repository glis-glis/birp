//
// Created by caduffm on 4/13/23.
//

#ifndef BIRP_TDETERMINISTIC_H
#define BIRP_TDETERMINISTIC_H

#include "TDistributions.h"
#include "TLikelihood.h"

//------------------------------------------
// TDeterministic
//------------------------------------------

class TDeterministic {
	// red box
protected:
	// gamma and times of change
	stattools::TParameterTyped<TypeGamma, 1> *_gamma        = nullptr;
	stattools::TParameterTyped<TypeTime, 1> *_timesOfChange = nullptr;

	TDistributionBase *_distr;

	// number of epochs
	size_t _numEpochs = 0;

	// rho
	coretools::TMatrix<size_t> _rho;

	// Jeffrey's prior
	double _oldJP = 0.0;

	void _fillRho(const TUniqueContainer<size_t> &Timepoints);
	std::vector<size_t> _getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const;
	double _calculateLogPhi_k(size_t k) const;

public:
	TDeterministic(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	               stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	virtual ~TDeterministic() = default;

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints);

	std::vector<double> calculatePhi(size_t j) const;

	void updateGamma(TLikelihood &LikelihoodCalculator);
};

#endif // BIRP_TDETERMINISTIC_H
