//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"
#include "TData.h"
#include "TEpochModel.h"

using namespace coretools::instances;

TBirpPrior::TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                       stattools::TObservationTyped<TypeEffort, 2> *Effort, const std::string &Filename)
    : _gamma(Gamma), _effort(Effort), _filename(Filename) {
	addPriorParameter({_gamma, _effort});

	// calculate nu
	_nu.resize(_numMethods);
	const auto &counts = *_storageBelow[0];
	for (size_t i = 0; i < _numMethods; ++i) {
		_nu[i] = 0.0;
		for (size_t k = 0; k < _numTimePoints; k++) { _nu[i] += (counts)(i, k); }
	}
}

std::string TBirpPrior::name() const { return "birp"; }

void TBirpPrior::initializeInferred() {
	// read times of change and fill rho
	// TEpochModel::_readTimesOfChange();
	const auto counts = _storageBelow[0];
	_numMethods       = counts->dimensions()[0];
	_numTimePoints    = counts->dimensions()[1];

	// initialize gamma accordingly
	_gamma->initStorage({_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});
}

void TBirpPrior::estimateInitialPriorParameters() {
	// TODO: implement

	// intitialize oldLL
	//_oldLL           = _calculateLL();
	//_oldLogJeffPrior = _calculateLogJeffreysPrior();
}

double TBirpPrior::_calculateLL(const coretools::TMatrix<double> &p) const {
	const auto &counts = *_storageBelow[0];

	double LL = 0.0;
	for (size_t k = 0; k < _numTimePoints; k++) {
		for (size_t i = 0; i < _numMethods; i++) { LL += (double)(counts)(i, k) * log(p(i, k)); }
	}

	return LL;
}

double TBirpPrior::_calculateLL() const {
	coretools::TMatrix<double> p; // = _calculateP();
	_calculateLL(p);
}

double TBirpPrior::_calculateLogJeffreysPrior(const coretools::TMatrix<double> &p) const {
	// TODO: implement

	return 0.0;
}

void TBirpPrior::_evaluateUpdateGamma(size_t m) {
	// calculate probabilities p for new gamma
	coretools::TMatrix<double> p; // = _calculateP();

	const double newLL   = _calculateLL(p);
	const double LLRatio = newLL - _oldLL;

	const double newLogJeffPrior = _calculateLogJeffreysPrior(p);
	const double logPriorRatio   = newLogJeffPrior - _oldLogJeffPrior;

	const double logH = LLRatio + logPriorRatio;
	if (_gamma->acceptOrReject(m, logH)) {
		// accepted -> update LL and jeffrey's Prior
		_oldLL           = newLL;
		_oldLogJeffPrior = newLogJeffPrior;
	}
}

void TBirpPrior::updateParams() {
	for (size_t m = 0; m < _numEpochs; m++) {
		_gamma->updateSpecificIndex(m);
		_evaluateUpdateGamma(m);
	}
}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {
	// LL is defined on all gamma jointly -> just return LL
	return _calculateLL(); // or return oldLL?
}

void TBirpPrior::_simulateUnderPrior(Storage *Counts) {
	// todo: implement

	// write as proper file formats
	TData::write(_filename, *Counts, _effort->storage());
}

double TBirpPrior::rho(size_t Time, size_t Epoch) const { return _rho(Time, Epoch); }
coretools::TMatrix<size_t> TBirpPrior::rhoMat() const { return _rho; }

bool TBirpPrior::_checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index, const std::string &Title) {
	return Infile.headerAt(index) == Title;
}

void TBirpPrior::_checkHeader(const coretools::TInputFile &Infile) {
	std::vector<std::string> headerVec{"species", "location", "method", "time", "counts", "efforts"};
	for (size_t i = 0; i < headerVec.size(); i++) {
		if (!_checkHeaderPlacement(Infile, i, headerVec[i])) {
			UERROR("Error in header of file", Infile.name(), " , please make sure header at column ", i + 1, " is ",
			       headerVec[i], ". Found: '", Infile.headerAt(i), "'");
		}
	}
}
void TBirpPrior::_readFile(const coretools::TInputFile &Infile) {
	
}

///
