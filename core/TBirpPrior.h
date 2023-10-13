//
// Created by caduffm on 7/20/22.
//

#ifndef TBIRPPRIOR_H
#define TBIRPPRIOR_H

#include "BirpTypes.h"
#include "TBirpConfiguration.h"
#include "TData.h"
#include "TModels.h"
#include "coretools/Files/TFile.h"
#include "coretools/Files/TOutputFile.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"

class TBirpPrior : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2> {
private:
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::Storage;
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::UpdatedStorage;

protected:
	// data
	TData &_data;

	// model
	std::vector<std::unique_ptr<TModelBase>> _models;

	// locations, timepoints, covNames: these are references to the objects that are stored in TBirpCore
	const TUniqueContainer<std::string> &_locations;
	const TUniqueContainer<TypeTime> &_timepoints;
	const TUniqueContainer<std::string> &_covariateEffortNames;
	const TUniqueContainer<std::string> &_covariateDetectionNames;

	size_t _numSpecies = 0;

	// prefix for writing
	std::string _out;

	// simulate values
	std::vector<double> _readNBarFromCommandLine() const;
	void _simulateUnderPrior(Storage *Data) override;

public:
	TBirpPrior() = default;
	TBirpPrior(TData &Data, const TUniqueContainer<std::string> &Locations,
	           const TUniqueContainer<TypeTime> &TimePoints,
	           const TUniqueContainer<std::string> &CovariateEffortNames,
	           const TUniqueContainer<std::string> &CovariateDetectionNames, size_t numSpecies,
	           TBirpCommonModel &CommonModel, std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn,
	           std::unique_ptr<TBirpStochasticAddOn> &StochasticAddOn, std::string Prefix);
	~TBirpPrior() override = default;

	void initializeInferred() override;

	[[nodiscard]] std::string name() const override;

	// estimate initial values
	void estimateInitialPriorParameters() override;

	// full log densities
	double getSumLogPriorDensity(const Storage &) const override;

	// update all parameters
	void updateParams() override;

	// write summaries of gamma posterior
	void burninHasFinished() override;
	void writeSummaryGammaPosterior(const std::string &Prefix);

	// getters for locations, methods, timepoints, species for tests
	const TUniqueContainer<std::string> &locations() const;
	const TUniqueContainer<TypeTime> &timepoints() const;
	const TData &data() const;

	stattools::TUpdateBase *getPtrToJointUpdater() { return _models[0]->getPtrToJointUpdater(); }
};

#endif // TBIRPPRIOR_H
