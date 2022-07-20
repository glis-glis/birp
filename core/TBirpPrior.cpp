//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"

TBirpPrior::TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                       stattools::TParameterTyped<TypeEffort, 1> *Effort) {}

std::string TBirpPrior::name() const { return "birp"; }

void TBirpPrior::initializeStorageOfPriorParameters() {}

void TBirpPrior::estimateInitialPriorParameters() {}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {}

void TBirpPrior::updateParams() {}

void TBirpPrior::_simulateUnderPrior(Storage *Data) {}
