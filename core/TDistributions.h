//
// Created by caduffm on 4/13/23.
//

#ifndef BIRP_TDISTRIBUTIONS_H
#define BIRP_TDISTRIBUTIONS_H

#include "BirpTypes.h"
#include "TData.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>

//---------------------------------------
// TDistribution
//---------------------------------------
class TDistributionBase {
protected:
	size_t _speciesIDinUniqueContainer = 0; // refers to ID in unique container //TODO add to constructor

public:
	TDistributionBase();
	virtual ~TDistributionBase() = default;

	virtual void initializeInferred(const TUniqueContainer<size_t> &Timepoints,
	                                const TUniqueContainer<std::string> &Locations, const TData &Data) = 0;
	virtual void update(const TData &Data)                                                             = 0;
	virtual void estimateInitialPriorParameters(const TData &Data)                                     = 0;

	virtual double calculateLL(const std::vector<double> &sigma, const TLocations &loc, const std::vector<double> &p,
	                           size_t index_in_counts, size_t ix_linNumMethLoc) const = 0;

	virtual void simulate(TData &methods, size_t N0,
	                      std::tuple<std::vector<std::string>, std::vector<std::string>,
	                                 std::vector<coretools::Probability>> &detectionProbs,
	                      const TUniqueContainer<std::string> &locations, bool SimulateCovariates) = 0;
};

//---------------------------------------
// TPoissonDistribution
//---------------------------------------
class TPoissonDistribution : public TDistributionBase {
protected:
public:
	TPoissonDistribution();

	~TPoissonDistribution() override = default;

	double calculateLL(const std::vector<double> &sigma, const TLocations &loc, const std::vector<double> &p,
	                   size_t index_in_counts, size_t ix_linNumMethLoc) const override;

	// simulation
	void simulate(TData &methods, size_t N0,
	              std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<coretools::Probability>>
	                  &detectionProbs,
	              const TUniqueContainer<std::string> &locations, bool SimulateCovariates) override;
};

//--------------------------------
// TNegBinDistribution
//--------------------------------

class TNegBinDistribution : public TDistributionBase {
protected:
	// parameters: mu and a
	stattools::TParameterTyped<TypeMu, 1> *_mu     = nullptr;
	stattools::TParameterTyped<TypeLogA, 1> *_logA = nullptr;

	// update functions
	void _updateMu(const TData &Data);
	void _updateLogA(const TData &Data);

	// initialization
	void _estimateInitialMu(const TData &Data);

	// simulate
	double _setAForSimulation(const TData &Data);

public:
	TNegBinDistribution(stattools::TParameterTyped<TypeMu, 1> *Mu, stattools::TParameterTyped<TypeLogA, 1> *LogA);

	~TNegBinDistribution() override = default;

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints, const TUniqueContainer<std::string> &Locations,
	                        const TData &Data) override;
	virtual void update(const TData &Data);
	void estimateInitialPriorParameters(const TData &Data) override;

	double calculateLL(const std::vector<double> &sigma, const TLocations &loc, const std::vector<double> &p,
	                   size_t index_in_counts, size_t ix_linNumMethLoc) const override;

	// simulate
	void simulate(TData &Data, size_t N0,
	              std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<coretools::Probability>>
	                  &detectionProbs,
	              const TUniqueContainer<std::string> &locations, bool SimulateCovariates) override;
};

#endif // BIRP_TDISTRIBUTIONS_H
