//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"
#include "TData.h"

TBirpPrior::TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                       stattools::TObservationTyped<TypeEffort, 2> *Effort, const std::string &Filename)
    : _gamma(Gamma), _effort(Effort), _filename(Filename) {
	addPriorParameter({_gamma, _effort});
}

std::string TBirpPrior::name() const { return "birp"; }

void TBirpPrior::_fillRho(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times) {
	_numEpochs = TimesOfChange.size() + 1;
	_rhos.zeros(Times.size(), _numEpochs);
	_rhos(0, 0) = 1;
}

void TBirpPrior::initializeStorageOfPriorParameters() {}

void TBirpPrior::estimateInitialPriorParameters() {}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const { return 0.0; }

void TBirpPrior::updateParams() {}

void TBirpPrior::_simulateUnderPrior(Storage *Counts) {
	// todo: implement

	// write as proper file formats
	TData::write(_filename, *Counts, _effort->storage());
}
