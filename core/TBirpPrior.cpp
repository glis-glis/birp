//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"
#include "TEpochModel.h"

using namespace coretools::instances;

TBirpPrior::TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma) { addPriorParameter(Gamma); }

std::string TBirpPrior::name() const { return "birp"; }

void TBirpPrior::initializeInferred() {
	// TODO: implement
	std::string filenameCounts = parameters().getParameterFilename("counts");
	_readFile(filenameCounts);

	// _epochModel.initializeGamma();
}

void TBirpPrior::estimateInitialPriorParameters() {
	// TODO: implement
}

void TBirpPrior::updateParams() {
	// TODO: implement
	// _epochModel.updateGamma();
}

double TBirpPrior::getSumLogPriorDensity(const Storage &) const {
	// LL is defined on all gamma jointly -> just return LL
	// TODO: implement
	// return _epochModel.calculateLL(); // or return oldLL?
	return 0.0;
}

void TBirpPrior::_simulateUnderPrior(Storage *) {
	// todo: implement
	// TODO: write as proper file formats
}

//
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
void TBirpPrior::_readFile(const std::string &Filename) {
	coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header);
	_checkHeader(infile);

	std::vector<std::string> vec;
	while (infile.read(vec)) {
		size_t locationIndex  = _locations.add(vec[1]);
		size_t methodIndex    = _methods.add(vec[2]);
		size_t timepointIndex = _timepoints.add(coretools::str::convertStringCheck<size_t>(vec[3]));
		TDataPoint pt(timepointIndex, coretools::str::convertStringCheck<size_t>(vec[4]),
		              coretools::str::convertStringCheck<double>(vec[5]));
		// loop over all species and check if it already exists
		for (size_t sp = 0; sp < _species.size(); sp++) {
			if (_species[sp].name() == vec[0]) { // did find it
				_species[sp].add(pt, locationIndex, methodIndex);
			} else { // did not find it
				_species.emplace_back(vec[0], pt, locationIndex, methodIndex);
			}
		}
	}
}
