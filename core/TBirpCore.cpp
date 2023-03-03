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
    return {"species", "location", "method", "time", "counts"};
}

void TBirpCore::_checkHeader(const coretools::TInputFile &Infile, bool hasCounts) {
    std::vector<std::string> headerVec = _getHeader();
    if (!hasCounts) { headerVec.pop_back(); }
    // std::cout << "headerVec: " << headerVec[0] << " " << headerVec[1] << " " << headerVec[2] << " " << headerVec[3]
    // << " " << headerVec[4] << std::endl;
    for (size_t i = 0; i < headerVec.size(); i++) {
        if (!_checkHeaderPlacement(Infile, i, headerVec[i])) {
            UERROR("Error in header of file", Infile.name(), " , please make sure header at column ", i + 1, " is ",
                   headerVec[i], ". Found: '", Infile.headerAt(i), "'");
        }
    }
}

void TBirpCore::_readForInference() {
    std::string filenameCounts = parameters().getParameterFilename("counts");
    _readCountFile(filenameCounts, true);
    _sortTimepoints();

    // TODO -> read effort
}

void TBirpCore::_readCountFile(const std::string &Filename, bool hasCounts) {
    std::string delim = "\t";
    std::string_view ending = coretools::str::readAfterLast(Filename, '.');
    if (ending == "csv") { delim = ","; }
    coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header, delim);
    if (hasCounts && infile.numCols() != 5) {
        UERROR("Expected file ", infile.name(), " to have 5 columns, but found ", infile.numCols(), ".");
    } else if (!hasCounts && infile.numCols() != 4) {
        UERROR("Expected file ", infile.name(), " to have 4 columns, but found ", infile.numCols(), ".");
    }
    _checkHeader(infile, hasCounts);

    std::vector<std::string> vec;
    while (infile.read(vec)) {
        size_t locationIndex = _locations.add(vec[1]);
        size_t methodIndex = _methods.add(vec[2]);
        size_t timepointIndex = _timepoints.add(coretools::str::fromStringCheck<size_t>(vec[3]));
        size_t counts = 0;
        if (hasCounts) { counts = coretools::str::fromStringCheck<size_t>(vec[4]); }
        TDataPoint pt(timepointIndex, counts, 0);  //TODO write new constructor
        // loop over all species and check if it already exists
        bool speciesExists = false;
        for (auto &species: _species) {
            if (species.name() == vec[0]) { // did find it
                speciesExists = true;
                species.add(pt, locationIndex, methodIndex);
                break;
            }
        }
        if (!speciesExists) { _species.emplace_back(vec[0], pt, locationIndex, methodIndex); }
    }
}

void TBirpCore::_readEffortFile(const std::string &Filename) {
    // 1) open file, check header
    std::string delim = "\t";
    std::string_view ending = coretools::str::readAfterLast(Filename, '.');
    if (ending == "csv") { delim = ","; }
    coretools::TInputFile infile(Filename, coretools::TFile_Filetype::header, delim);
    // TODO how to check header? check that first three cols are correct + at least 4 cols.
    size_t numCov = infile.numCols() - 3;  // the first three columns are not covariates

    // 2) fill _covariateNames from corresponding header column
    for (size_t i = 0; i < numCov; ++i) {
        _covariateNames.add(infile.headerAt(i + 3));
    }

    // 3) for each row, check if this combination of meth/loc/timepoint exists in TSpecies, if no: go to next row
    std::vector<std::string> row;
    while (infile.read(row)) {
        size_t locationIndex = _locations.getIndex(row[0]);
        size_t methodIndex = _methods.getIndex(row[1]);
        size_t timepointIndex = _timepoints.getIndex(coretools::str::fromStringCheck<size_t>(row[2]));
        bool exists = false;
        for (auto &species: _species) {
            if (species.IJK_exists(methodIndex, locationIndex, timepointIndex)) {
                exists = true;
                break;
            };
        }
        if (!exists) { continue; } // will go to next row in file
        // store in TCovariates
        for (size_t i = 0; i < numCov; ++i) {
            TCovariate_icjk cov(timepointIndex,
                                coretools::str::fromStringCheck<TypeCovariate>(row[i + 3])); //converts to TypeCovariate
            _covariates.addTimePoint(cov, methodIndex, i, locationIndex);
        }

    }
    // TODO 6) loop over all species and meth/locs & timepoints: check if that combination of ijk exists in TCovariates, throw error if not
    for (auto &species: _species){
        species.checkCovExists();
    }
    // TODO: timepoints inside vector will not be sorted, but will have the correct timepoint_id -> is this a problem?




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
    auto header = _getHeader();
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
