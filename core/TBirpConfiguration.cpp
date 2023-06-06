//
// Created by caduffm on 4/18/23.
//

#include "TBirpConfiguration.h"

#include "TBirpPrior.h"
#include "TData.h"
#include "coretools/Types/strongTypes.h"
#include "coretools/devtools.h"
#include "stattools/MCMC/TMCMC.h"
#include "stattools/Priors/TPriorDirichlet.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorNormal.h"

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

//------------------------------------------
// TBirpCommonModel
//------------------------------------------

void TBirpCommonModel::_createAlpha(size_t i, TData &Data, const TUniqueContainer<std::string> &CovariateEffortNames,
                                    const std::string &Filename) {
	// create alpha
	auto name  = "alpha_" + Data[i].name();
	auto prior = std::make_shared<stattools::prior::TDirichletFixed<stattools::TParameterBase, TypeAlpha, 1>>();
	TParameterDefinition def(Filename, readAlphaDirichlet("dirichlet_alpha", CovariateEffortNames.size(),
	                                                      10.0)); // alpha = 10 is quite arbitrary choice
	alpha[i] = std::make_shared<stattools::TParameterTyped<TypeAlpha, 1>>(name, prior, def);
}

void TBirpCommonModel::_createBeta0(size_t i, TData &Data, const std::string &Filename) {
	// create beta
	auto name  = "beta0_" + Data[i].name();
	auto prior = std::make_shared<stattools::prior::TNormalFixed<stattools::TParameterBase, TypeBeta, 1>>();
	TParameterDefinition def(Filename, "0,1");
	beta0[i] = std::make_shared<stattools::TParameterTyped<TypeBeta, 1>>(name, prior, def);
}

void TBirpCommonModel::_createBeta(size_t i, TData &Data, const std::string &Filename) {
	// create beta
	auto name  = "beta_" + Data[i].name();
	auto prior = std::make_shared<stattools::prior::TNormalFixed<stattools::TParameterBase, TypeBeta, 1>>();
	TParameterDefinition def(Filename, "0,1");
	beta[i] = std::make_shared<stattools::TParameterTyped<TypeBeta, 1>>(name, prior, def);
}

TBirpCommonModel::TBirpCommonModel(TData &Data, const TUniqueContainer<std::string> &CovariateEffortNames,
                                   const std::string &Filename, stattools::TDAGBuilder &DAGBuilder)
    // TODO: change to improper prior!
    : gamma("gamma", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeGamma, 1>>(),
            {Filename}),
      timesOfChange("timesOfChange",
                    std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeTime, 1>>(),
                    {Filename}) {

	// create alpha and beta
	const size_t numMethods = Data.size();
	alpha.resize(numMethods);
	beta.resize(numMethods);
	beta0.resize(numMethods);
	for (size_t i = 0; i < numMethods; ++i) {
		_createAlpha(i, Data, CovariateEffortNames, Filename);
		_createBeta0(i, Data, Filename);
		_createBeta(i, Data, Filename);
	}

	gamma.getDefinition().setJumpSizeForAll(false);

	// add to DAG
	DAGBuilder.addToDAG({&gamma, &timesOfChange});
	for (auto &a : alpha) { DAGBuilder.addToDAG(a.get()); }
	for (auto &b : beta0) { DAGBuilder.addToDAG(b.get()); }
	for (auto &b : beta) { DAGBuilder.addToDAG(b.get()); }
}

//------------------------------------------
// TBirpNegBinomAddOn
//------------------------------------------

TBirpNegBinomAddOn::TBirpNegBinomAddOn(TData &Data, const std::string &Filename, stattools::TDAGBuilder &DAGBuilder)
    : sigmaMu("sigmaMu",
              std::make_shared<stattools::prior::TExponentialFixed<stattools::TParameterBase, TypeSigmaMu, 1>>(),
              {Filename, "1"}),
      mu("mu",
         std::make_shared<stattools::prior::TDirichletVarInferred<stattools::TParameterBase, TypeMu, 1, TypeSigmaMu>>(
             &sigmaMu),
         {Filename, readAlphaDirichlet("alpha_mu", Data.numMethLoc(0),
                                       1.0)}), // todo speciesID 0 because right now single species
      logA("logA", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogA, 1>>(),
           {Filename}) {

	// add to DAG
	DAGBuilder.addToDAG({&mu, &sigmaMu, &logA});
}

//------------------------------------------
// TBirpStochasticAddOn
//------------------------------------------

TBirpStochasticAddOn::TBirpStochasticAddOn(const std::string &Filename, stattools::TDAGBuilder &DAGBuilder)
    : logPhi("logPhi", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogPhi, 2>>(),
             {Filename}),
      logSigma("logSigma",
               std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogSigma, 1>>(),
               {Filename}) {

	logPhi.getDefinition().setJumpSizeForAll(false);
	logPhi.getDefinition().setEqualNumberOfUpdates(false); // first logPhi is not updated

	TypeLogSigma ::setMin(-10.); // prevent numeric issues with very small sigma: 1/sigma can result in nan

	// add to DAG
	DAGBuilder.addToDAG({&logPhi, &logSigma});
}

//------------------------------------------
// TBirpModel
//------------------------------------------

TBirpModel::TBirpModel(TData &Data, const TUniqueContainer<std::string> &Locations,
                       const TUniqueContainer<std::size_t> &TimePoints,
                       const TUniqueContainer<std::string> &SpeciesNames,
                       const TUniqueContainer<std::string> &CovariateEffortNames,
                       const TUniqueContainer<std::string> &CovariateDetectionNames, const std::string &Filename,
                       stattools::TDAGBuilder &DAGBuilder)
    : _commonModel(Data, CovariateEffortNames, Filename, DAGBuilder) {

	if (parameters().parameterExists("negativeBinomial") && parameters().parameterExists("stochastic")) {
		_negBinomAddOn   = std::make_unique<TBirpNegBinomAddOn>(Data, Filename, DAGBuilder);
		_stochasticAddOn = std::make_unique<TBirpStochasticAddOn>(Filename, DAGBuilder);
	} else if (parameters().parameterExists("negativeBinomial")) {
		_negBinomAddOn = std::make_unique<TBirpNegBinomAddOn>(Data, Filename, DAGBuilder);
	} else if (parameters().parameterExists("stochastic")) {
		_stochasticAddOn = std::make_unique<TBirpStochasticAddOn>(Filename, DAGBuilder);
	}

	// create birp prior (Poisson or Negative Binomial)
	_birpPrior =
	    std::make_shared<TBirpPrior>(Data, Locations, TimePoints, CovariateEffortNames, CovariateDetectionNames,
	                                 SpeciesNames.size(), _commonModel, _negBinomAddOn, _stochasticAddOn, Filename);

	// create fake storage: data are stored in TSpecies, but stattools requires TMultiDimensionalStorage
	coretools::TMultiDimensionalStorage<TypeCounts, 2> fakeData;
	_data = std::make_unique<stattools::TObservationTyped<TypeCounts, 2>>("likelihood", _birpPrior, fakeData,
	                                                                      stattools::TObservationDefinition(Filename));

	// add to DAG
	DAGBuilder.addToDAG(_data.get());
	DAGBuilder.buildDAG();
}

void TBirpModel::writeSummaryGammaPosterior(const std::string &Prefix) {
	_birpPrior->writeSummaryGammaPosterior(Prefix);
}
