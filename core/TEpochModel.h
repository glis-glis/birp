//
// Created by liam on 23/08/22.
//

#ifndef BIRP_TEPOCHMODEL_H
#define BIRP_TEPOCHMODEL_H

#include "BirpTypes.h"
#include "TMatrix.h"
#include "TObservationTyped.h"
#include "TParameterTyped.h"

class TEpochModel public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2> {
private:
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::Storage;
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::UpdatedStorage;

	// Birp parameters and co-variates
	stattools::TParameterTyped<TypeGamma, 1> *_gamma     = nullptr;
	stattools::TObservationTyped<TypeEffort, 2> *_effort = nullptr;

};

#endif // BIRP_TEPOCHMODEL_H
