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
	                           size_t index_in_counts, size_t ix_linNumMethLoc, size_t ix_method) const = 0;

	virtual void simulate(TData &Data, const TModelBase &Model, const std::vector<double> &n_bar,
	                      bool SimulateCovariates) = 0;
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
	                   size_t index_in_counts, size_t ix_linNumMethLoc, size_t ix_method) const override;

	// simulation
	void simulate(TData &Data, const TModelBase &Model, const std::vector<double> &n_bar,
	              bool SimulateCovariates) override;
};

//--------------------------------
// TNegBinDistribution
//--------------------------------

class TNegBinDistribution : public TDistributionBase {
protected:
	// parameters: mu (or N) and b
	stattools::TParameterTyped<TypeMuOrN, 1> *_muOrN = nullptr;
	stattools::TParameterTyped<TypeB, 1> *_b         = nullptr;
	bool _isMu;

	double _getMuOrN(size_t j, size_t ij) const;

	// update functions
	void _updateMu(const TData &Data, TModelBase &Model);
	void _updateB(const TData &Data, TModelBase &Model);
	std::pair<double, std::vector<double>> _calcHastingsRatioN(const std::vector<size_t> &ij, const TData &Data,
	                                                           TModelBase &Model);
	void _updateN(const TData &Data, TModelBase &Model);

	// initialization
	void _initializeN(const TUniqueContainer<std::string> &Locations);
	void _initializeMu(const TData &Data, const TUniqueContainer<std::string> &Locations);
	void _estimateInitialMu(const TData &Data);
	void _estimateInitialN(const TData &Data);
	std::vector<double> _getNuPerIJ(const TData &Data);

	// simulate
	static std::vector<double> _readAFromCommandLine(const TData &Data);

public:
	TNegBinDistribution(size_t SpeciesIDinUniqueContainer, stattools::TParameterTyped<TypeMuOrN, 1> *MuOrN,
	                    stattools::TParameterTyped<TypeB, 1> *B);

	~TNegBinDistribution() override = default;

	void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations) override;
	void update(const TData &Data, TModelBase &Model) override;
	void estimateInitialParameters(const TData &Data) override;

	double calculateLL(const std::vector<double> &sigma, const TLocations &loc, const std::vector<double> &p,
	                   size_t index_in_counts, size_t ix_linNumMethLoc, size_t ix_method) const override;

	// simulate
	void simulate(TData &Data, const TModelBase &Model, const std::vector<double> &n_bar,
	              bool SimulateCovariates) override;
};

#endif // BIRP_TDISTRIBUTIONS_H
