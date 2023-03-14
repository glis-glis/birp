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

TBirpModel::TBirpModel(std::vector<TSpecies> *Species, TUniqueContainer<std::string> *Locations,
                       TUniqueContainer<std::string> *Methods, TUniqueContainer<std::size_t> *TimePoints,
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

    // currently only allow for single species
    if (Species->size() != 1) { UERROR("Currently only supports single species!"); }

    // create birp prior (Poisson or Negative Binomial)
    std::shared_ptr<TBirpPrior> birp;
    if (NegativeBinomial) {
        birp = std::make_shared<TBirpPrior>(Species, Locations, Methods, TimePoints, &_gamma, &_timesOfChange, &_mu,
                                            &_logA, Filename);
    } else {
        birp =
                std::make_shared<TBirpPrior>(Species, Locations, Methods, TimePoints, &_gamma, &_timesOfChange,
                                             Filename);
    }

    // create fake storage: data are stored in TSpecies, but stattools requires TMultiDimensionalStorage
    coretools::TMultiDimensionalStorage<TypeCounts, 2> fakeData;
    _data = std::make_unique<stattools::TObservationTyped<TypeCounts, 2>>("data", birp, fakeData,
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

void TBirpCore::_checkHeader(const coretools::TInputFile &Infile, bool hasCounts, std::vector<size_t> countIndices,
                             std::vector<size_t> effortIndices) {
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
                   " , should contain at least 4 columns: 'method', 'location' and at least two more columns starting with 'counts' and 'cov', respectively.");
        }
    } else {
        if (2 + countIndices.size() + effortIndices.size() < 3) {
            UERROR("Header of file", Infile.name(),
                   " , should contain at least 3 columns: 'method', 'location' and at least one column starting with 'cov'.");
        }
    }
}

void TBirpCore::_readForInference(bool hasCounts) {
    std::string filename = parameters().getParameterFilename("counts");
    _readFile(filename, hasCounts);
    _sortTimepoints();
}

std::vector<size_t> TBirpCore::getColumnIndex(coretools::TInputFile infile, std::string string) {
    std::vector<size_t> countIndices;
    for (size_t i = 0; i < infile.numCols(); i++) {
        if (coretools::str::stringStartsWith(infile.headerAt(i), string)) {
            countIndices.push_back(i);
        }
    }
    return countIndices;
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

    _checkHeader(infile, bool hasCounts, countIndices, effortIndices);

    std::vector<std::string> row;
    while (infile.read(row)) {
        size_t locationIndex = _locations.add(row[0]);
        size_t timepointIndex = _timepoints.add(coretools::str::fromStringCheck<size_t>(row[1]));

        // boot one vector together (stiflä) of all species counts, and one of all covariates
        std::vector<TypeCounts> countVec(countIndices.size());
        std::vector<TypeCovariate> covVec(effortIndices.size());
        for (size_t c = 0; c < countVec.size(); c++) {
            countVec[c] = (TypeCounts) row[countIndices[c]];
        }
        for (size_t c = 0; c < covVec.size(); c++) {
            countVec[c] = (TypeCovariate) row[effortIndices[c]];
        }
        // create timepoint
        TTimepoints pt(countVec, covVec);
        // add timepoint to data, either in existing object or create new one
        _methods[0].add(pt, locationIndex, timepointIndex); // TODO how to know which  method?

        // TODO ignoriere hasCounts for a start
    }
}


void TBirpCore::_setEffort_singleCovariate() {
    // if effort is fixed: write the corresponding values (for each location & timepoint) to TDataPoint
    for (auto &method: _covariates) {
        if (method.covariateIsEffort()) {
            auto cov = method.begin(); // we know there is only one, so no need for loop
            for (auto &loc: *cov) {
                for (auto &time: loc) {
                    for (auto &species: _species) {
                        // _species[method_id][location_id][time_id].effort = value; --> does this work?
                        species[][time.timepoint_id].effort = time.covariate;
                    }
                }
            }

        }
    } else {
        // do nothing for the moment
    }
}

void TBirpCore::_sortTimepoints() {
    // sort timepoints and corresponding dataseries for every species
    // check if data has more than one entry with effort > 0
    for (auto &species: _species) { species.checkData(_locations, _methods); }
    std::vector<size_t> sortingIndex = coretools::rankSort(_timepoints);
    _timepoints = coretools::sortContainerByRank(_timepoints, sortingIndex);
    for (auto &species: _species) { species.sorttimes(sortingIndex); }
}

void TBirpCore::infer() {
    // read data
    _readForInference();

    // create model
    stattools::TDAGBuilder dagBuilder;
    bool negBinom = parameters().parameterExists("negativeBinomial");
    TBirpModel model(&_species, &_locations, &_methods, &_timepoints, _outname, dagBuilder, negBinom);

    // run MCMC
    stattools::TMCMC mcmc;
    mcmc.runMCMC(_outname, &dagBuilder);
}

//-----------------------------------
// Simulate

void TBirpCore::_writeFile(const std::string &Filename) {
    auto header = _getHeaderCounts();
    coretools::TOutputFile outfile(Filename, header);
    // write file line by line
    for (auto &species: _species) { species.write(outfile, _locations, _methods, _timepoints); }
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

void TBirpCore::_fillMethodsFromCommandLine() {
    size_t numMethods = parameters().getParameterWithDefault("numMethods", 2);
    for (size_t m = 0; m < numMethods; ++m) { _methods.add("Method_" + coretools::str::toString(m)); }
}

void TBirpCore::simulate() {
    // read incomplete counts file (if provided)
    if (parameters().parameterExists("simulationFile")) {
        std::string filename = parameters().getParameterFilename("simulationFile");
        _readCountFile(filename, false);
    } else {
        // else read main dimensions from command line
        _fillTimepointsFromCommandLine();
        _fillLocationsFromCommandLine();
        _fillMethodsFromCommandLine();
        _species.emplace_back("Seal", _locations.size(), _methods.size(), _timepoints.size());
    }

    // create model
    stattools::TDAGBuilder dagBuilder;
    bool negBinom = parameters().parameterExists("negativeBinomial");
    TBirpModel model(&_species, &_locations, &_methods, &_timepoints, _outname, dagBuilder, negBinom);

    // run simulator
    stattools::TSimulator simulator;
    simulator.simulate(_outname, &dagBuilder);

    // write Birp file
    _writeFile(_outname + "_simulatedCounts.txt");
}
