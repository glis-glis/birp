//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"

TBirpPrior::TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                       stattools::TObservationTyped<TypeEffort, 2> *Effort)
    : _gamma(Gamma), _effort(Effort) {
	addPriorParameter({_gamma, _effort});
}

std::string TBirpPrior::name() const { return "birp"; }

void TBirpPrior::initializeStorageOfPriorParameters() {}

void TBirpPrior::estimateInitialPriorParameters() {}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {}

void TBirpPrior::updateParams() {}

void TBirpPrior::_simulateUnderPrior(Storage *Data) {}
