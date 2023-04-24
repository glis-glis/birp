//
// Created by caduffm on 4/18/23.
//

#ifndef BIRP_TBIRPCONFIGURATION_H
#define BIRP_TBIRPCONFIGURATION_H

#include "BirpTypes.h"
#include "TData.h"
#include "coretools/Main/TParameters.h"
#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorUniform.h"

//------------------------------------------
// TBirpCommonModel
//------------------------------------------

class TBirpCommonModel {
public:
	// all parameters of model (including prior parameters)
	stattools::TParameterTyped<TypeGamma, 1> gamma;
	stattools::TParameterTyped<TypeTime, 1> timesOfChange;

	// covariates: regression coefficients
	std::vector<std::shared_ptr<stattools::TParameterTyped<TypeSigmaBeta, 1>>> sigmaBeta;
	std::vector<std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>>> beta;

	TBirpCommonModel(TData &Data, const TUniqueContainer<std::string> &CovariateNames, const std::string &Filename,
	                 stattools::TDAGBuilder &DAGBuilder);
};

//------------------------------------------
// TBirpNegBinomAddOn
//------------------------------------------

class TBirpNegBinomAddOn {
public:
	// only relevant for negative binomial model
	stattools::TParameterTyped<TypeSigmaMu, 1> sigmaMu;
	stattools::TParameterTyped<TypeMu, 1> mu;
	stattools::TParameterTyped<TypeLogA, 1> logA;

	TBirpNegBinomAddOn(TData &Data, const std::string &Filename, stattools::TDAGBuilder &DAGBuilder);
};

//------------------------------------------
// TBirpStochasticAddOn
//------------------------------------------

class TBirpStochasticAddOn {
public:
	// only relevant for stochastic model
	stattools::TParameterTyped<TypeLogPhi, 2> logPhi;
	stattools::TParameterTyped<TypeSigma, 1> sigma;

	TBirpStochasticAddOn(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder);
};

//------------------------------------------
// TBirpModel
//------------------------------------------

class TBirpModel {
private:
	TBirpCommonModel _commonModel;
	std::unique_ptr<TBirpNegBinomAddOn> _negBinomAddOn;
	std::unique_ptr<TBirpStochasticAddOn> _stochasticAddOn;

	std::unique_ptr<stattools::TObservationTyped<TypeCounts, 2>> _data;

public:
	TBirpModel(TData &Data, const TUniqueContainer<std::string> &Locations,
	           const TUniqueContainer<std::size_t> &TimePoints, const TUniqueContainer<std::string> &SpeciesNames,
	           const TUniqueContainer<std::string> &CovariateNames, const std::string &Filename,
	           stattools::TDAGBuilder &DAGBuilder);
};

#endif // BIRP_TBIRPCONFIGURATION_H
