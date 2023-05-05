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

// forward declarations
class TModelBase;

//---------------------------------------
// TDistribution
//---------------------------------------
class TDistributionBase {
protected:
	size_t _speciesIDinUniqueContainer = 0; // refers to ID in unique container

public:
	TDistributionBase(size_t SpeciesIDinUniqueContainer);
	virtual ~TDistributionBase() = default;

	virtual void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations) = 0;
	virtual void update(const TData &Data, TModelBase &Model)                                  = 0;
	virtual void estimateInitialParameters(const TData &Data)                                  = 0;

	virtual double calculateLL(const std::vector<double> &sigma, const TLocations &loc, const std::vector<double> &p,
	                           size_t index_in_counts, size_t ix_linNumMethLoc) const = 0;

	virtual void simulate(TData &Data, const TModelBase &Model, size_t N0, bool SimulateCovariates) = 0;
};

//---------------------------------------
// TPoissonDistribution
//---------------------------------------
class TPoissonDistribution : public TDistributionBase {
public:
	TPoissonDistribution(size_t SpeciesIDinUniqueContainer);
	~TPoissonDistribution() override = default;

	void initialize(const TData &, const TUniqueContainer<std::string> &) override{};
	void update(const TData &, TModelBase &) override{};
	void estimateInitialParameters(const TData &) override{};

	double calculateLL(const std::vector<double> &sigma, const TLocations &loc, const std::vector<double> &p,
	                   size_t index_in_counts, size_t ix_linNumMethLoc) const override;

	// simulation
	void simulate(TData &Data, const TModelBase &Model, size_t N0, bool SimulateCovariates) override;
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
	void _updateMu(const TData &Data, TModelBase &Model);
	void _updateLogA(const TData &Data, TModelBase &Model);

	// initialization
	void _estimateInitialMu(const TData &Data);

	// simulate
	double _setAForSimulation(const TData &Data);

public:
	TNegBinDistribution(size_t SpeciesIDinUniqueContainer, stattools::TParameterTyped<TypeMu, 1> *Mu,
	                    stattools::TParameterTyped<TypeLogA, 1> *LogA);

	~TNegBinDistribution() override = default;

	void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations) override;
	void update(const TData &Data, TModelBase &Model) override;
	void estimateInitialParameters(const TData &Data) override;

	double calculateLL(const std::vector<double> &sigma, const TLocations &loc, const std::vector<double> &p,
	                   size_t index_in_counts, size_t ix_linNumMethLoc) const override;

	// simulate
	void simulate(TData &Data, const TModelBase &Model, size_t N0, bool SimulateCovariates) override;
};

#endif // BIRP_TDISTRIBUTIONS_H
