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

class TDeterministic : public TModelBase {
	// red box
protected:
	// gamma and times of change
	TGamma *_gamma = nullptr;

	// distribution: Poisson or negative binomial
	TDistributionBase *_distr = nullptr;

	// store temporary phi
	std::vector<TypePhi> _phi;

	void _fillPhi();
	double _calculateLogPhi_k(size_t k) const;

public:
	TDeterministic(TGamma *Gamma);
	~TDeterministic() override = default;

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints);

	TypePhi phi(size_t location_id, size_t timepoint_id) const override;

	void updateGamma(const TData &Data, TLikelihood &Likelihood);
};

#endif // BIRP_TDETERMINISTIC_H
