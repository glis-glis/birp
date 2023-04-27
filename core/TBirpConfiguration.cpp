//
// Created by caduffm on 4/18/23.
//

#include "TBirpConfiguration.h"

#include "TBirpCore.h"
#include "TBirpPrior.h"
#include "TData.h"
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

//------------------------------------------
// TBirpCommonModel
//------------------------------------------

TBirpCommonModel::TBirpCommonModel(TData &Data, const TUniqueContainer<std::string> &CovariateNames,
                                   const std::string &Filename, stattools::TDAGBuilder &DAGBuilder)
    // TODO: change to improper prior!
    : gamma("gamma", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeGamma, 1>>(),
            {Filename}),
      timesOfChange("timesOfChange",
                    std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeTime, 1>>(),
                    {Filename}) {

	// create beta
	sigmaBeta.resize(Data.size());
	beta.resize(Data.size());
	for (size_t i = 0; i < Data.size(); ++i) {
		// create sigma beta (variance of Dirichlet prior on beta_i)
		auto name_sigma_beta = "sigma_beta_" + Data[i].name();
		auto prior_sigma_beta =
		    std::make_shared<stattools::prior::TExponentialFixed<stattools::TParameterBase, TypeSigmaBeta, 1>>();
		TParameterDefinition def_sigma_beta({Filename, "1"});
		if (!Data[i].inferEffort()) { def_sigma_beta.update(false); }
		sigmaBeta[i] = std::make_shared<stattools::TParameterTyped<TypeSigmaBeta, 1>>(name_sigma_beta, prior_sigma_beta,
		                                                                              def_sigma_beta);

		// create beta
		auto name  = "beta_" + Data[i].name();
		auto prior = std::make_shared<
		    stattools::prior::TDirichletVarInferred<stattools::TParameterBase, TypeBeta, 1, TypeSigmaBeta>>(
		    sigmaBeta[i].get());
		TParameterDefinition def(Filename, readAlphaDirichlet("alpha_beta", CovariateNames.size(), 1.0));
		beta[i] = std::make_shared<stattools::TParameterTyped<TypeBeta, 1>>(name, prior, def);
	}

	gamma.getDefinition().setJumpSizeForAll(false);

	// add to DAG
	DAGBuilder.addToDAG({&gamma, &timesOfChange});
	for (auto &b : beta) { DAGBuilder.addToDAG(b.get()); }
	for (auto &s : sigmaBeta) { DAGBuilder.addToDAG(s.get()); }
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
      sigma("sigma", std::make_shared<stattools::prior::TExponentialFixed<stattools::TParameterBase, TypeSigma, 1>>(),
            {Filename, "10"}) {

	// add to DAG
	DAGBuilder.addToDAG({&logPhi, &sigma});
}

//------------------------------------------
// TBirpModel
//------------------------------------------

TBirpModel::TBirpModel(TData &Data, const TUniqueContainer<std::string> &Locations,
                       const TUniqueContainer<std::size_t> &TimePoints,
                       const TUniqueContainer<std::string> &SpeciesNames,
                       const TUniqueContainer<std::string> &CovariateNames, const std::string &Filename,
                       stattools::TDAGBuilder &DAGBuilder)
    : _commonModel(Data, CovariateNames, Filename, DAGBuilder) {

	if (parameters().parameterExists("negativeBinomial") && parameters().parameterExists("stochastic")) {
		_negBinomAddOn   = std::make_unique<TBirpNegBinomAddOn>(Data, Filename, DAGBuilder);
		_stochasticAddOn = std::make_unique<TBirpStochasticAddOn>(Filename, DAGBuilder);
	} else if (parameters().parameterExists("negativeBinomial")) {
		_negBinomAddOn = std::make_unique<TBirpNegBinomAddOn>(Data, Filename, DAGBuilder);
	} else if (parameters().parameterExists("stochastic")) {
		_stochasticAddOn = std::make_unique<TBirpStochasticAddOn>(Filename, DAGBuilder);
	}

	// create birp prior (Poisson or Negative Binomial)
	std::shared_ptr<TBirpPrior> birp =
	    std::make_shared<TBirpPrior>(Data, Locations, TimePoints, CovariateNames, SpeciesNames.size(), _commonModel,
	                                 _negBinomAddOn, _stochasticAddOn, Filename);

	// create fake storage: data are stored in TSpecies, but stattools requires TMultiDimensionalStorage
	coretools::TMultiDimensionalStorage<TypeCounts, 2> fakeData;
	_data = std::make_unique<stattools::TObservationTyped<TypeCounts, 2>>("fakeData", birp, fakeData,
	                                                                      stattools::TObservationDefinition(Filename));

	// add to DAG
	DAGBuilder.addToDAG(_data.get());
	DAGBuilder.buildDAG();
}
