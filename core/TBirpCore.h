/*
 * TBirpCore.h
 *
 *  Created on: Jul 20, 2022
 *      Author: mcaduff
 */

#ifndef TEXAMPLETASK_H_
#define TEXAMPLETASK_H_

#include "BirpTypes.h"
#include "TBirpPrior.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TTask.h"
#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorUniform.h"
#include "TData.h"

//------------------------------------------
// TBirpModel
//------------------------------------------

class TBirpModel {
private:
    // all parameters of model (including prior parameters)
    stattools::TParameterTyped<TypeGamma, 1> _gamma;
    stattools::TParameterTyped<TypeTime, 1> _timesOfChange;

    // only relevant for negative binomial model
    stattools::TParameterTyped<TypeSigmaMu, 1> _sigmaMu;
    stattools::TParameterTyped<TypeMu, 1> _mu;
    stattools::TParameterTyped<TypeLogA, 1> _logA;

    std::unique_ptr<stattools::TObservationTyped<TypeCounts, 2>> _data;

public:
    TBirpModel(std::vector<TMethods> *Methods, TUniqueContainer<std::string> *Locations,
               TUniqueContainer<std::size_t> *TimePoints, TUniqueContainer<std::string> *speciesNames,
               const std::string &Filename, stattools::TDAGBuilder &DAGBuilder, bool NegativeBinomial);
};

//--------------------------------------
// TBirpCore
//--------------------------------------

class TBirpCore {
protected:
    std::string _outname;

    // locations, methods, timepoints
    TUniqueContainer<std::string> _locations;
    std::vector<TMethods> _methods;
    TUniqueContainer<size_t> _timepoints;
    TUniqueContainer<std::string> _covariateNames;
    TUniqueContainer<std::string> _speciesNames;

    // read input files
    [[nodiscard]] std::vector<std::string> _getHeader() const;

    [[nodiscard]] bool _checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index, const std::string &Title) const;

    void _checkHeader(const coretools::TInputFile &Infile, bool hasCounts, const std::vector<size_t> &countIndices,
                      const std::vector<size_t> &effortIndices);

    std::vector<size_t>
    getColumnIndex(const coretools::TInputFile &infile, const std::string &string);
    std::vector<std::string> _getSpeciesNames(const coretools::TInputFile &infile, const std::vector<size_t>& countIndices);
    void _readFile(const std::string &Filename, bool hasCounts);

    void _sortTimepoints();

    void _readForInference(bool hasCounts);

    // simulate
    std::vector<size_t> _readN0() const;

    void _fillTimepointsFromCommandLine();

    void _fillLocationsFromCommandLine();


    // write
    std::vector<std::string> _getMethodspecificHeader(const TMethods& method);
    void _writeFile();

public:
    TBirpCore();

    // main tasks
    void infer();

    void simulate();
};

//--------------------------------------
// Tasks
//--------------------------------------

class TTask_infer : public coretools::TTask {
public:
    TTask_infer() { _explanation = "Inferring rates of change in species densities from count data."; };

    void run() override {
        TBirpCore birp;
        birp.infer();
    };
};

class TTask_simulate : public coretools::TTask {
public:
    TTask_simulate() { _explanation = "Simulating count data under Birp model."; };

    void run() override {
        TBirpCore birp;
        birp.simulate();
    };
};

#endif /* TEXAMPLETASK_H_ */
