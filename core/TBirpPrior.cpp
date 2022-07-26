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
}

std::string TBirpPrior::name() const { return "birp"; }

std::vector<double> TBirpPrior::_getEpochStarts(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times) {
    double maxTimeOfChange;
    maxTimeOfChange = max_element(TimesOfChange.begin(), TimesOfChange.end());
    double maxTime = std::max(Times[Times.size()], maxTimeOfChange + 1);
    //
    std::vector<double> _epochStarts = 0;
    _epochStarts.insert(_epochStarts.end(), TimesOfChange.begin(), TimesOfChange.end());
    _epochStarts.push_back(maxTime);
    return _epochStarts;
}

 void TBirpPrior::initializeStorageOfPriorParameters() {
	// read times of change and fill rho
	_readTimesOfChange();
	const auto counts = _storageBelow[0];
	_numMethods       = counts->dimensions()[0];
	_numTimePoints    = counts->dimensions()[1];

	// initialize gamma accordingly
	_gamma->initializeStorageBasedOnPrior({_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});
}

coretools::TMatrix<double> TBirpPrior::_fillRho(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times) {
	_numEpochs = TimesOfChange.size() + 1;
	_rho.zeros((int)Times.size(), (int)_numEpochs);
    std::vector<double> epochStartT = TBirpPrior::_getEpochStarts(TimesOfChange, Times);
	// fill rhos
	_rho(0, 0) = 1.0;
    for(auto e = 1; e > _numEpochs; ++e){
        _rho(Times >= epochStartT[e+1],e) = epochStartT[e+1] - epochStartT[e];
        bool withinEpoch = Times < epochStartT[e+1] & Times > epochStartT[e];
        _rho(withinEpoch,e) = Times[withinEpoch] - epochStartT[e];
    }
    return _rho;
}

void TBirpPrior::initializeStorageOfPriorParameters() {}


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

void TBirpPrior::estimateInitialPriorParameters() {
	// TODO: implement

	// intitialize oldLL
	_oldLL           = _calculateLL(_storageBelow[0]);
	_oldLogJeffPrior = _calculateLogJeffreysPrior();
}

double TBirpPrior::_calculatePhi_k(size_t k) const {
	// TODO: implement eq (1) (for a specific time point k)
	return 1.0;
}

std::vector<double> TBirpPrior::_calculatePhi() const {
	std::vector<double> phi(_numTimePoints);
	for (size_t k = 0; k < _numTimePoints; k++) { phi[k] = _calculatePhi_k(k); }
	return phi;
}

double TBirpPrior::_calculateSigma_ik(size_t i, size_t k, double Phi_k) const {
	return Phi_k * _effort->value(_effort->getIndex({i, k}));
}

coretools::TMatrix<double> TBirpPrior::_calculateP() const {
	std::vector<double> phi = _calculatePhi();

	coretools::TMatrix<double> p((int)_numMethods, (int)_numTimePoints);
	for (size_t i = 0; i < _numMethods; i++) {
		double sum = 0.0;
		for (size_t k = 0; k < _numTimePoints; k++) {
			// calculate sigma and fill into p (will normalize later)
			const double sigma = _calculateSigma_ik(i, k, phi[k]);
			p(i, k)            = sigma;
			sum += sigma;
		}
		// now normalize by row sum!
		for (size_t k = 0; k < _numTimePoints; k++) { p(i, k) /= sum; }
	}
	return p;
}

double TBirpPrior::_calculateLL(Storage *Counts) const {
	// calculate probabilities p for new gamma
	coretools::TMatrix<double> p = _calculateP();

	double LL = 0.0;
	for (size_t i = 0; i < _numMethods; i++) {
		for (size_t k = 0; k < _numTimePoints; k++) { LL += (*Counts)(i, k) * log(p(i, k)); }
	}

	return LL;
}

double TBirpPrior::_calculateLogJeffreysPrior() const {
	// TODO: implement
	return 0.0;
}

void TBirpPrior::_evaluateUpdateGamma(size_t m, Storage *Counts) {
	const double newLL   = _calculateLL(Counts);
	const double LLRatio = newLL - _oldLL;

	const double newLogJeffPrior = _calculateLogJeffreysPrior();
	const double logPriorRatio   = newLogJeffPrior - _oldLogJeffPrior;

	const double logH = LLRatio + logPriorRatio;
	if (_gamma->acceptOrReject(m, logH)) {
		// accepted -> update LL and jeffrey's Prior
		_oldLL           = newLL;
		_oldLogJeffPrior = newLogJeffPrior;
	}
}

void TBirpPrior::updateParams() {
	auto counts = this->_storageBelow[0];

	for (size_t m = 0; m < _numEpochs; m++) {
		_gamma->updateSpecificIndex(m);
		_evaluateUpdateGamma(m, counts);
	}
}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {
	// LL is defined on all gamma jointly -> just return LL
	return _calculateLL(_storageBelow[0]);
}

void TBirpPrior::_simulateUnderPrior(Storage *Counts) {
	// todo: implement

	// write as proper file formats
	TData::write(_filename, *Counts, _effort->storage());
}

double TBirpPrior::rho(size_t Time, size_t Epoch) { return _rho(Time, Epoch); }
