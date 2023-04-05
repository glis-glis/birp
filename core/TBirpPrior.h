//
// Created by caduffm on 7/20/22.
//

#ifndef TBIRPPRIOR_H
#define TBIRPPRIOR_H

#include "BirpTypes.h"
#include "TData.h"
#include "TEpochModel.h"
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
	std::vector<TMethods> *_methods = nullptr;

	// model
	std::vector<std::unique_ptr<TEpochModel>> _epochModel;

	// locations, timepoints, covNames: these are pointers to the objects that are stored in TBirpCore
	TUniqueContainer<std::string> *_locations      = nullptr;
	TUniqueContainer<size_t> *_timepoints          = nullptr;
	TUniqueContainer<std::string> *_covariateNames = nullptr;

	size_t _numSpecies = 0;

	// prefix for writing
	std::string _out;

	// simulate values
	std::vector<size_t> _readN0FromCommandLine() const;
	coretools::TMultiDimensionalStorage<double, 2> _readDetectionProbFile() const;
	void _simulateUnderPrior(Storage *Data) override;
	void _simulateFromFile();
	void _simulateFromCommandLine();

public:
	TBirpPrior() = default;

	TBirpPrior(std::vector<TMethods> *Methods, TUniqueContainer<std::string> *Locations,
	           TUniqueContainer<std::size_t> *TimePoints, TUniqueContainer<std::string> *CovariateNames,
	           size_t numSpecies, stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	           stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
	           std::vector<std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>>> &Beta, std::string Prefix);

	TBirpPrior(std::vector<TMethods> *Methods, TUniqueContainer<std::string> *Locations,
	           TUniqueContainer<std::size_t> *TimePoints, TUniqueContainer<std::string> *CovariateNames,
	           size_t numSpecies, stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	           stattools::TParameterTyped<TypeTime, 1> *TimesOfChange, stattools::TParameterTyped<TypeMu, 1> *Mu,
	           stattools::TParameterTyped<TypeLogA, 1> *LogA,
	           std::vector<std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>>> &Beta, std::string Prefix);

	~TBirpPrior() override = default;

	void initializeInferred() override;

	[[nodiscard]] std::string name() const override;

	// estimate initial values
	void estimateInitialPriorParameters() override;

	// full log densities
	double getSumLogPriorDensity(const Storage &) const override;

	// update all parameters
	void updateParams() override;

	// getters for locations, methods, timepoints, species for tests
	const TUniqueContainer<std::string> &locations() const;

	const TUniqueContainer<size_t> &timepoints() const;

	const std::vector<TMethods> &methods() const;
};

#endif // TBIRPPRIOR_H
