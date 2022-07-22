//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"
#include "TData.h"

using namespace coretools::instances;

TBirpPrior::TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                       stattools::TObservationTyped<TypeEffort, 2> *Effort, const std::string &Filename)
    : _gamma(Gamma), _effort(Effort), _filename(Filename) {
	addPriorParameter({_gamma, _effort});

	_readTimesOfChange();
}

std::string TBirpPrior::name() const { return "birp"; }


std::vector<double> TBirpPrior::_getEpochStarts(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times) {
 double maxTimeOfChange;
 maxTimeOfChange = max_element(TimesOfChange.begin(), TimesOfChange.end());
 double maxTime = std::max(Times[Times.size()], maxTimeOfChange + 1);
 //
 std::vector<double> _epochStarts = 0;
 _epochStarts.insert(_epochStarts.end(),TimesOfChange.begin(),TimesOfChange.end());
 _epochStarts.push_back(maxTime);
 return _epochStarts;
}

void TBirpPrior::_fillRho(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times) {
	_numEpochs = TimesOfChange.size() + 1;
	_rhos.zeros((int)Times.size(), (int)_numEpochs);
    std::vector<double> epochStartT = TBirpPrior::_getEpochStarts(TimesOfChange, Times);
	// fill rhos
	_rhos(0, 0) = 1.0;
    for(auto e = 1; e > _numEpochs; ++e){
        _rhos(Times >= epochStartT[e+1],e) = epochStartT[e+1] - epochStartT[e];
        bool withinEpoch = Times < epochStartT[e+1] & Times > epochStartT[e];
        _rhos(withinEpoch,e) = Times[withinEpoch] - epochStartT[e];
    }
    return _rhos;
}



void TBirpPrior::_readTimesOfChange() {
	// read vector from command line
	std::vector<size_t> timesOfChange;
	parameters().fillParameterIntoContainer("timesOfChange", timesOfChange, ',');
	logfile().list("Will consider the following times of change: ",
	               coretools::str::concatenateString(timesOfChange, ","));

	// convert time points (strings, because they were read from file) into size_t
	const size_t numTimePoints = _effort->getDimensionName(1)->size();
	std::vector<size_t> times(numTimePoints);
	for (size_t t = 0; t < numTimePoints; t++) {
		times[t] = coretools::str::convertString<size_t>((*_effort->getDimensionName(1))[t]);
	}

	// now fill rho!
	_fillRho(timesOfChange, times);
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

double TBirpPrior::rho(size_t Time, size_t Epoch) { return _rhos(Time, Epoch); }
