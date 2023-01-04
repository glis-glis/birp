//
// Created by caduffm on 7/20/22.
//

#include "TBirpPrior.h"
#include "TEpochModel.h"
#include "coretools/Storage/TDataFile.h"

using namespace coretools::instances;

TBirpPrior::TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma) {
    addPriorParameter(Gamma);
    _out = parameters().getParameterWithDefault("out", "birp");
}

std::string TBirpPrior::name() const { return "birp"; }

void TBirpPrior::initializeInferred() {
    // TODO: implement
    read();
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

std::vector<std::string> TBirpPrior::_getHeader() const {
    return {"species", "location", "method", "time", "effort", "counts"};
}

void TBirpPrior::_checkHeader(const coretools::TInputFile &Infile) {
    std::vector<std::string> headerVec = _getHeader();
    for (size_t i = 0; i < headerVec.size(); i++) {
        if (!_checkHeaderPlacement(Infile, i, headerVec[i])) {
            UERROR("Error in header of file", Infile.name(), " , please make sure header at column ", i + 1, " is ",
                   headerVec[i], ". Found: '", Infile.headerAt(i), "'");
        }
    }
}

void TBirpPrior::read() {
    std::string filenameCounts = parameters().getParameterFilename("counts");
    _readFile(filenameCounts, true);
    _sortTimepoints();
}

void TBirpPrior::_readFile(const std::string &Filename, bool hasCounts) {
    std::string delim = "\t";
    std::string_view ending = coretools::str::readAfterLast(Filename, '.');
    if (ending == "csv") { delim = ","; }
    coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header, delim);
    // TODO if it has counts: 6 columns, otherwise 5.
    _checkHeader(infile);

    std::vector<std::string> vec;
    while (infile.read(vec)) {
        size_t locationIndex = _locations.add(vec[1]);
        size_t methodIndex = _methods.add(vec[2]);
        size_t timepointIndex = _timepoints.add(coretools::str::fromStringCheck<size_t>(vec[3]));
        size_t counts = 0;
        if (hasCounts) { counts = coretools::str::fromStringCheck<size_t>(vec[5]); }
        TDataPoint pt(timepointIndex, counts, coretools::str::fromStringCheck<double>(vec[4]));
        // loop over all species and check if it already exists
        bool speciesExists = false;
        for (size_t sp = 0; sp < _species.size(); sp++) {
            if (_species[sp].name() == vec[0]) { // did find it
                speciesExists = true;
                _species[sp].add(pt, locationIndex, methodIndex);
                break;
            }
        }
        if (!speciesExists) { _species.emplace_back(vec[0], pt, locationIndex, methodIndex); }
    }
}

// sort timepoints and corresponding dataseries for every species
void TBirpPrior::_sortTimepoints() {
    // check if data has more than one entry with effort > 0
    for (size_t sp = 0; sp < _species.size(); sp++) { _species[sp].checkData(locations(), methods()); }
    std::vector<size_t> sortingIndex = coretools::rankSort(_timepoints);
    _timepoints = coretools::sortContainerByRank(_timepoints, sortingIndex);
    for (size_t sp = 0; sp < _species.size(); sp++) { _species[sp].sorttimes(sortingIndex); }
}

void TBirpPrior::_writeFile(const std::string &Filename) {
    auto header = _getHeader();
    coretools::TOutputFile outfile(Filename, header);
    // write file line by line
    for (auto &species: _species) {
        for (auto &lm: species) {
            // lm: loop over all location-methods combinations
            std::string location = _locations[lm.location_id()];
            std::string method = _methods[lm.method_id()];
            for (auto &t: lm) { // loop over all time points
                size_t timepoint = _timepoints[t.timepoint_id];
                outfile << species.name() << location << method << timepoint << t.effort << t.counts;
                outfile.endln();
            }
        }
    }
}

void TBirpPrior::_simulateFromFile() {
    auto detectionProbs = _readDetectionProbFile();
    auto N0 = _readN0();
    std::string filename = parameters().getParameterFilename("simulationFile");
    _readFile(filename, false);
    for (size_t i = 0; i < _epochModel.size(); i++) {
        _epochModel[i].simulate(_species[i], N0[i], detectionProbs, _locations, _methods);
    }
    _writeFile(_out + "_simulated.txt");
}

void TBirpPrior::_simulateFromCommandLine() {
    auto timepoints= _readTimepoints();
	// TODO how to fill a single number?
    std::vector<size_t> numLoc;
    std::vector<size_t> numLocDefault = {2};
    parameters().fillParameterIntoContainerWithDefault("numLocations", numLoc, ',',numLocDefault);
    std::vector<size_t> numMeth;
    std::vector<size_t> numMethDefault = {2};
    parameters().fillParameterIntoContainerWithDefault("numMethods", numMeth, ',',numMethDefault);
    size_t numLocations = numLoc[0];
    size_t numMethods = numMeth[0];
    _species.emplace_back("Moose", numLocations, numMethods, timepoints.size());
	for (size_t i = 0; i < _epochModel.size(); i++) {
		_epochModel[i].simulate(_species[i]);
	}
	_writeFile(_out + "_simulated.txt");
}


// TODO are detection probabilities the same across all species or one file per species?
coretools::TMultiDimensionalStorage<double, 2> TBirpPrior::_readDetectionProbFile() {
    std::string filename = parameters().getParameterFilename("detectionProbFile");

    auto locationNames = std::make_shared<coretools::TNamesStrings>();
    auto methodNames = std::make_shared<coretools::TNamesStrings>();

    // read file
    coretools::TDataReaderBase<double, 2> reader;
    size_t guessNumRows = 50; // doesn't matter, not too many anyways
    coretools::TMultiDimensionalStorage<double, 2> detectionProbs;
    reader.read(filename, &detectionProbs, {locationNames, methodNames}, {'\n', '\t'}, {true, true}, {0, 0},
                {false, false}, guessNumRows, coretools::colNames_multiline);
    return detectionProbs;
    // TODO make default p=1
}

std::vector<size_t> TBirpPrior::_readN0() {
    std::vector<size_t> N0;
    std::vector<size_t> N0default = {100};
    parameters().fillParameterIntoContainerWithDefault("N0", N0, ',', N0default);
    return N0;
}

std::vector<size_t> TBirpPrior::_readTimepoints() {
    std::vector<size_t> timepoints;
    std::vector<size_t> timepointsDdefault = {1,2,3};
    parameters().fillParameterIntoContainerWithDefault("timepoints", timepoints, ',', timepointsDdefault);
    return timepoints;
}





// getters for test functions
const TUniqueContainer<std::string> &TBirpPrior::locations() const { return _locations; }

const TUniqueContainer<std::string> &TBirpPrior::methods() const { return _methods; }

const TUniqueContainer<size_t> &TBirpPrior::timepoints() const { return _timepoints; }

const std::vector<TSpecies> &TBirpPrior::species() const { return _species; }

///
