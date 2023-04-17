//
// Created by caduffm on 4/13/23.
//

#ifndef BIRP_TDETERMINISTIC_H
#define BIRP_TDETERMINISTIC_H

#include "TDistributions.h"
#include "TGamma.h"
#include "TLikelihood.h"

//------------------------------------------
// TDeterministic
//------------------------------------------

class TDeterministic {
	// red box
protected:
	// gamma and times of change
	TGamma *_gamma = nullptr;

	TDistributionBase *_distr = nullptr;

	double _calculateLogPhi_k(size_t k) const;

public:
	TDeterministic(TGamma *Gamma);
	virtual ~TDeterministic() = default;

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints);

	std::vector<double> calculatePhi(size_t j) const;

	void updateGamma(TLikelihood &LikelihoodCalculator);
};

#endif // BIRP_TDETERMINISTIC_H
