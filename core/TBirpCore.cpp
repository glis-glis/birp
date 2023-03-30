/*
 * TBirpCore.cpp
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#include "TBirpCore.h"
#include "coretools/Types/strongTypes.h"
#include "stattools/MCMC/TMCMC.h"
#include "stattools/Priors/TPriorDirichlet.h"
#include "stattools/Priors/TPriorExponential.h"
#include "TData.h"

using namespace coretools::instances;
using namespace stattools;
using namespace stattools::prior;

std::string readAlphaDirichlet(const std::string &ArgName, size_t Size, coretools::StrictlyPositive<double> Default) {
    // read values of alpha from command line
    std::vector<coretools::StrictlyPositive<double>> alphas;
    if (parameters().parameterExists(ArgName)) {
        parameters().fillParameterIntoContainer(ArgName, alphas, ',');
    } else {
        alphas.push_back(Default);
    }

    if (alphas.empty()) { UERROR("Invalid argument '", ArgName, "': is empty!"); }
    if (alphas.size() == 1) { return coretools::str::repeat(alphas[0], Size); }
    if (alphas.size() == Size) { return coretools::str::concatenateString(alphas, ","); }

    UERROR("Size of argument '", ArgName, "' (", alphas.size(), ") does not match expected size (", Size, ").");
}

TBirpModel::TBirpModel(std::vector<TMethods> *Methods, TUniqueContainer<std::string> *Locations,
                       TUniqueContainer<std::size_t> *TimePoints, TUniqueContainer<std::string> *speciesNames,
                       const std::string &Filename, stattools::TDAGBuilder &DAGBuilder, bool NegativeBinomial)
// TODO: change to improper prior!
        : _gamma("gamma", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeGamma, 1>>(),
                 {Filename}),
          _timesOfChange("timesOfChange",
                         std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeTime, 1>>(),
                         {Filename}),
          _sigmaMu("sigmaMu",
                   std::make_shared<stattools::prior::TExponentialFixed<stattools::TParameterBase, TypeSigmaMu, 1>>(),
                   {Filename, "1"}),
          _mu("mu",
              std::make_shared<stattools::prior::TDirichletVarInferred<stattools::TParameterBase, TypeMu, 2, TypeSigmaMu>>(
                      &_sigmaMu),
              {Filename, readAlphaDirichlet("alpha_mu", Methods->size() * Locations->size(), 1.0)}),
          _logA("logA", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogA, 1>>(),
                {Filename}) {


    // create birp prior (Poisson or Negative Binomial)
    std::shared_ptr<TBirpPrior> birp;
    if (NegativeBinomial) {
        birp = std::make_shared<TBirpPrior>(Methods, Locations, TimePoints, speciesNames, &_gamma, &_timesOfChange, &_mu,
                                            &_logA, Filename);
    } else {
        birp =
                std::make_shared<TBirpPrior>(Methods, Locations, TimePoints, speciesNames, &_gamma, &_timesOfChange,
                                             Filename);
    }

    // create fake storage: data are stored in TSpecies, but stattools requires TMultiDimensionalStorage
    coretools::TMultiDimensionalStorage<TypeCounts, 2> fakeData;
    _data = std::make_unique<stattools::TObservationTyped<TypeCounts, 2>>
            ("fakeData", birp, fakeData,
             stattools::TObservationDefinition(Filename));

    // add to DAG
    DAGBuilder.addToDAG({&_gamma, &_timesOfChange});
    if (NegativeBinomial) { DAGBuilder.addToDAG({&_mu, &_sigmaMu, &_logA}); }
    DAGBuilder.addToDAG(_data.get());
    DAGBuilder.buildDAG();
};

//--------------------------------------
// TBirpCore
//--------------------------------------

TBirpCore::TBirpCore() {
    using namespace coretools::instances;
    if (parameters().parameterExists("out")) {
        _outname = parameters().getParameter<std::string>("out");
        logfile().list("Will write files with prefix '", _outname, "'. (parameter 'out')");
    } else {
        _outname = "birp";
        logfile().list("Will write files with default prefix '", _outname,
                       "'. (You can change this using the parameter 'out')");
    }
}


//-----------------------------------
// Infer

bool TBirpCore::_checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index,
                                      const std::string &Title) const {
    return Infile.headerAt(index) == Title;
}

std::vector<std::string> TBirpCore::_getHeader() const {
    return {"location", "timepoint"};
}

void
TBirpCore::_checkHeader(const coretools::TInputFile &Infile, bool hasCounts, const std::vector<size_t> &countIndices,
                        const std::vector<size_t> &effortIndices) {
    std::vector<std::string> headerVec = _getHeader();
    // check whether first two columns are correct
    for (size_t i = 0; i < headerVec.size(); i++) {
        if (!_checkHeaderPlacement(Infile, i, headerVec[i])) {
            UERROR("Error in header of file", Infile.name(), " , please make sure header at column ", i + 1, " is ",
                   headerVec[i], ". Found: '", Infile.headerAt(i), "'");
        }
    }
    // check if number of columns is at least 4 (with counts) or 3 (without counts)
    if (hasCounts) {
        if (2 + countIndices.size() + effortIndices.size() < 4) {
            UERROR("Header of file", Infile.name(),
                   " , should contain at least 4 columns: 'location', 'timepoint' and at least two more columns starting with 'counts_' and 'cov_', respectively.");
        }
    } else {
        if (2 + countIndices.size() + effortIndices.size() < 3) {
            UERROR("Header of file", Infile.name(),
                   " , should contain at least 3 columns: 'method', 'location' and at least one column starting with 'cov'.");
        }
    }
}

void TBirpCore::_readForInference(bool hasCounts) {
    std::string input = parameters().getParameterFilename("data");
    bool dataIsFile = std::filesystem::exists(
            input);  // if multiple files are supplied separated by commas, then this will be false
    std::vector<std::string> filenames;
    if (dataIsFile) {
        coretools::TInputFile infile(input, coretools::TFile_Filetype::header);
        // first check if maybe the file is the actual counts file and not a list of files
        if (infile.numCols() != 1) {
            coretools::str::fillContainerFromString(input, filenames, ',');
        } else { // if not, then its a list of files (methods), so we fill them in
            coretools::TInputFile infile(input, coretools::TFile_Filetype::fixed); // fixed: no header
            if (infile.numCols() != 1) {
                UERROR("Expected only one column in file ", input, " .");
            }
            std::vector<std::string> row;
            while (infile.read(row)) {
                filenames.push_back(row[0]);
            }
        }
    } else {
        coretools::str::fillContainerFromString(input, filenames, ',');
    }
    for (auto name: filenames) {
        _readFile(name, hasCounts);
    }
    _sortTimepoints();
}

std::vector<size_t> TBirpCore::getColumnIndex(const coretools::TInputFile &infile, const std::string &string) {
    std::vector<size_t> countIndices;
    for (size_t i = 0; i < infile.numCols(); i++) {
        if (coretools::str::stringStartsWith(infile.headerAt(i), string)) {
            countIndices.push_back(i);
        }
    }
    return countIndices;
}

std::vector<std::string>
TBirpCore::_getSpeciesNames(const coretools::TInputFile &infile, const std::vector<size_t> &countIndices) {
    std::vector<std::string> names(countIndices.size());
    std::string test = infile.headerAt(2);
    for (size_t s = 0; s < countIndices.size(); s++) {
        names[s] = coretools::str::readAfter(infile.headerAt(countIndices[s]), '_');
    }
    return names;
}

void TBirpCore::_readFile(const std::string &Filename, bool hasCounts) {
    std::string delim = "\t";
    std::string_view ending = coretools::str::readAfterLast(Filename, '.');
    if (ending == "csv") { delim = ","; }
    coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header, delim);

    // get position of columns that contain counts: might be an empty vector in case the file has no counts
    std::vector<size_t> countIndices = getColumnIndex(infile, "counts");

    // get position of columns that contain covariates (efforts)
    std::vector<size_t> effortIndices = getColumnIndex(infile, "cov");

    _checkHeader(infile, hasCounts, countIndices, effortIndices);

    std::vector<std::string> speciesNames = _getSpeciesNames(infile, countIndices);
    std::vector<size_t> species_id = {};
    for (auto name: speciesNames){
        species_id.push_back(_speciesNames.add(name));
    }

    // add new method: name is derived from filename
    std::string_view methodNameFull = coretools::str::readAfterLast(Filename, '/');
    std::string_view methodNameCropped = coretools::str::readAfter(methodNameFull, '_');
    std::string_view methodName = coretools::str::readBeforeLast(methodNameCropped, '.');
    _methods.emplace_back(methodName, _methods.size(), species_id); // the index of the method is the current size of the _methods vector

    std::vector<std::string> row;
    while (infile.read(row)) {
        size_t location_id = _locations.add(row[0]);
        size_t timepoint_id = _timepoints.add(coretools::str::fromStringCheck<size_t>(row[1]));

        // boot one vector together (stiflä) of all species counts, and one of all covariates
        std::vector<TypeCounts> countVec(countIndices.size());
        std::vector<TypeCovariate> covVec(effortIndices.size());
        for (size_t c = 0; c < countVec.size(); c++) {
            countVec[c] = (TypeCounts) row[countIndices[c]];
        }
        for (size_t c = 0; c < covVec.size(); c++) {
            covVec[c] = (TypeCovariate) row[effortIndices[c]];
        }
        // create timepoint
        TTimepoints pt(countVec, covVec, timepoint_id);
        // add timepoint to data, either in existing object or create new one
        _methods.back().add(pt, location_id); // back: goes to last element

        // TODO ignoriere hasCounts for a start
    }
    for (auto &m: _methods){
        m.fillIDVectors(_speciesNames.size());
    }
}


void TBirpCore::_sortTimepoints() {
    // sort timepoints and corresponding locations(?) for every method/species
    // check if data (location) has more than one entry with effort > 0

    std::vector<size_t> sortingIndex = coretools::rankSort(_timepoints);
    _timepoints = coretools::sortContainerByRank(_timepoints, sortingIndex);
    for (auto &method: _methods) {
        method.checkEffortBiggerThanZero();
        method.sorttimes(sortingIndex);
    }
}

void TBirpCore::infer() {
    // read data
    _readForInference(true);

    // create model
    stattools::TDAGBuilder dagBuilder;
    bool negBinom = parameters().parameterExists("negativeBinomial");
    TBirpModel model(&_methods, &_locations, &_timepoints, &_speciesNames, _outname, dagBuilder, negBinom);

    // run MCMC
    stattools::TMCMC mcmc;
    mcmc.runMCMC(_outname, &dagBuilder);
}

//-----------------------------------
// Simulate

std::vector<std::string> TBirpCore::_getMethodspecificHeader(const TMethods &method) {
    std::vector<std::string> header = _getHeader();
    // add species name to base header vector
    for (size_t i = 0; i < _speciesNames.size(); i++) {
        header.push_back("counts_" + _speciesNames[i]);
    }
    // add covariate name to base header vector
    for (size_t i = 0; i < method.covariateNames().size(); i++) {
        header.push_back("cov_" + method.covariateNames()[i]);
    }
    return header;
}

void TBirpCore::_writeFile() {
    // make second file that stores names of all methodFiles generated during the simulation
    coretools::TOutputFile fileList("allFilesGenerated");
    // write file line by line
    for (auto &method: _methods) {
        // get filename here.
        std::string Filename = _outname + "_" + method.methodName() + "_simulatedCounts.txt";
        // get header here
        std::vector<std::string> header = _getMethodspecificHeader(method);
        // append to fileList
        fileList << method.methodName() << coretools::endl; //TODO where is file?
        // now create actual file with data and write
        coretools::TOutputFile outfile(Filename, header);
        method.write(outfile, _locations, _timepoints);
    }
}

void TBirpCore::_fillTimepointsFromCommandLine() {
    std::vector<size_t> timepoints;
    std::vector<size_t> timepointsDefault = {1, 2, 3};
    parameters().fillParameterIntoContainerWithDefault("timepoints", timepoints, ',', timepointsDefault);

    // fill into unique container
    for (const auto &timepoint: timepoints) { _timepoints.add(timepoint); }
}

void TBirpCore::_fillLocationsFromCommandLine() {
    size_t numLocations = parameters().getParameterWithDefault("numLocations", 2);
    for (size_t l = 0; l < numLocations; ++l) { _locations.add("Location_" + coretools::str::toString(l)); }
}


void TBirpCore::simulate() {
    // read incomplete counts file (if provided)
    if (parameters().parameterExists("simulationFile")) {
        std::string filename = parameters().getParameterFilename("simulationFile");
        _readFile(filename, false);
    } else {
        // else read main dimensions from command line
        _fillTimepointsFromCommandLine();
        _fillLocationsFromCommandLine();
        size_t numMethods = parameters().getParameterWithDefault("numMethods", 1);
        size_t numSpecies = parameters().getParameterWithDefault("numSpecies", 1);
        size_t numCovariates = parameters().getParameterWithDefault("numCovariates", 1);
        // create vector of speciesnames species_1, species_2 ... as many as there are. currently only supports one species
        std::vector<size_t> species_ids(numSpecies);
        for (size_t n = 0; n < numSpecies; n++) {
            species_ids[n] = _speciesNames.add("species_" + coretools::str::toString(n + 1));
        }
        std::vector<std::string> covariateNames(numCovariates);
        for (size_t n = 0; n < numCovariates; n++) {
            covariateNames[n] = coretools::str::toString(n + 1);
        }
        for (size_t m = 0; m < numMethods; ++m) {
            _methods.emplace_back("Method_" + coretools::str::toString(m + 1), _locations.size(),
                                  _timepoints.size(), m, species_ids, covariateNames);
        }
        for (auto &m: _methods){
            m.fillIDVectors(_speciesNames.size());
        }
    }

    // create model
    stattools::TDAGBuilder dagBuilder;
    bool negBinom = parameters().parameterExists("negativeBinomial");
    TBirpModel model(&_methods, &_locations, &_timepoints, &_speciesNames, _outname, dagBuilder, negBinom);

    // run simulator
    stattools::TSimulator simulator;
    simulator.simulate(_outname, &dagBuilder); // goes to ... DAG ... TBirpPrior

    // write Birp file
    _writeFile();
}
